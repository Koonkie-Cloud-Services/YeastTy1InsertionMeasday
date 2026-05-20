#!/usr/bin/env python

import os
import argparse
import sys
import re
from typing import List, Dict
from collections import namedtuple

import gffutils
from gffutils import inspect
from gffutils import feature as feature_mod
from gffutils import create

FEATURE_TYPE_NAME_MAP = {
    "gene": "gene",
    "telomere": "telomere",
    "subtelomere": "subtelomere",
    "pericentromere": "pericentromere",
    "centromere": "centromere",
    "origin_of_replication": "origin_of_replication",
    "tRNA": "tRNA",
    "snRNA": "snRNA",
    "rRNA": "rRNA",
    "pseudogene": "pseudogene",
    "sequence_feature": "sequence_feature",  # Highly abundant
    "mobile_genetic_element": "mobile_element",  # Highly abundant
    "long_terminal_repeat": "LTR",
    "ncRNA": "other_RNA",
    "snoRNA": "snoRNA",
    "antisense_RNA": "other_RNA",
    "telomerase_RNA": "other_RNA",
    "SRP_RNA": "other_RNA",
    "RNase_P_RNA": "other_RNA",
    "RNase_MRP_RNA": "other_RNA",
    "regulatory_region": "regulatory_region",
}
SKIPPED_FEATURE_TYPES = [
    "mRNA",
    "CDS",
    "exon",
    "region",
    "transcript"
]

# The following feature types are replicated on both
# strands to ensure that they are covered by anti-features
# regardless of the strand of the original feature
STRAND_REPLICATED_FEATURE_TYPES = list(FEATURE_TYPE_NAME_MAP.keys())

# RNA feature types that take precedence over their parent gene entries
RNA_FEATURE_TYPES = {
    'tRNA', 'snoRNA', 'snRNA', 'rRNA', 'ncRNA',
    'antisense_RNA', 'telomerase_RNA', 'SRP_RNA',
    'RNase_P_RNA', 'RNase_MRP_RNA'
}

VALID_FEATURE_TYPES = set(FEATURE_TYPE_NAME_MAP)
VALID_FEATURE_TYPES.update(set(SKIPPED_FEATURE_TYPES))
TMP_GFF_PATH = "tmp.gff"

BUFFER_SPEC = namedtuple("BUFFER_SPEC", ["feature_name", "prefix", "sep", "length"])

MGE_RE_LIST = [
    re.compile(r"retrotransposon:(\w+)(Ty\d)-(\d+)$")
]

def get_options(args=sys.argv[1:]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert the features in a GFF file to BED format"
    )
    parser.add_argument(
        "--gff_in",
        dest="p_i_gff",
        type=str,
        help="The GFF file to convert"
    )
    parser.add_argument(
        "--annotation_table",
        dest="p_i_annot",
        type=str,
        required=False,
        help="A BED file containing additional annotations to add to the GFF"
    )
    parser.add_argument(
        "--bed_out",
        dest="p_o_bed",
        type=str,
        help="The name of the BED file to create"
    )
    parser.add_argument(
        "--buffer_features",
        dest="p_i_buffer_features",
        type=str,
        required=False,
        help="A TSV file containing features to add surrounding upstream and downstream region features."
        "The table should have the following columns: 1. Feature name, 2. Buffer prefix, 3. Separator character 4. Buffer length (bp)"
    )
    return parser.parse_args(args=args)


def create_gff_db(gff_file: str) -> gffutils.FeatureDB:
    # Inspect the GFF file
    res = inspect.inspect(gff_file)
    sys.stdout.write("\nFeature types in GFF file:\n")
    for ft, val in sorted(res["featuretype"].items(), key=lambda item: item[1], reverse=True):
        if ft not in VALID_FEATURE_TYPES:
            raise ValueError(f"Unknown feature type: {ft}")
        else:
            print(f"\t{ft}\t{val}")

    # Load the GFF file to a database
    gff_db = create.create_db(
        gff_file,
        dbfn=":memory:",
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
        force_merge_fields=[
            'seqid',
            'source',
            'featuretype'
        ]
    )
    return gff_db


def get_genome_of_feature(feature: gffutils.Feature) -> str:
    try:
        return feature.attributes["genome"][0]
    except KeyError:
        return None

def verify_feature_type(feature: gffutils.Feature) -> bool:
    return feature.featuretype in VALID_FEATURE_TYPES

def format_feature_as_bed(feature: gffutils.Feature) -> str:
    bed_row = [
        feature.chrom,
        feature.start,
        feature.stop,
        feature.id,
        feature.strand
    ]
    return "\t".join([str(i) for i in bed_row])

def should_replicate_on_both_strands(feature: gffutils.Feature) -> bool:
    """Check if a feature type should be replicated on both strands"""
    return feature.featuretype in STRAND_REPLICATED_FEATURE_TYPES

def create_opposite_strand_feature(feature: gffutils.Feature) -> gffutils.Feature:
    """Create a copy of the feature with the opposite strand
    
    The feature ID and all other properties remain the same, only the strand is flipped.
    """
    opposite_strand = '-' if feature.strand == '+' else '+'
    opposite_feature = feature_mod.feature_from_line(
        f"{feature.chrom}\t{feature.source}\t{feature.featuretype}\t{feature.start}\t{feature.stop}\t.\t{opposite_strand}\t.\tID={feature.id}"
    )
    opposite_feature.id = feature.id
    opposite_feature.featuretype = feature.featuretype
    return opposite_feature

def parse_gff_to_bed(gff_db: gffutils.FeatureDB, bed_file: str, **kwargs) -> None:
    outf = open(bed_file, "w")

    chrom_prefix = kwargs.get("chrom_prefix", "chr")
    genome_skip = kwargs.get("genome_skip", ['mitochondrion', 'plasmid'])
    chrom_skip = kwargs.get("chrom_skip", ['MT'])
    feature_types_to_skip = kwargs.get("feature_types_to_skip", SKIPPED_FEATURE_TYPES)
    chrom_skip.extend([chrom_prefix + chrom for chrom in chrom_skip])

    annot_features = kwargs.get("annot_features", {})
    if annot_features is None:
        annot_features = {}

    buffer_spec_index = kwargs.get("buffer_spec_index", {})

    # Identify genes that are parents of RNA features to avoid writing
    # duplicate entries (RNA annotations take precedence over gene entries)
    genes_with_rna_children = set()
    for rna_type in RNA_FEATURE_TYPES:
        for rna_feat in gff_db.all_features(featuretype=rna_type):
            for parent in gff_db.parents(rna_feat, featuretype='gene'):
                genes_with_rna_children.add(parent.id)

    chrom = ""
    acc = 0
    chrom_id_map = {}
    chrom_lengths = {}
    for ft_inst in gff_db.all_features():
        acc += 1

        if ft_inst.featuretype == "region":
            # Set the chromosome name for the following features
            try:
                chrom = ft_inst.attributes.get("Name")[0]
                chrom_id_map[ft_inst.chrom] = chrom_prefix + chrom
                chrom_lengths[chrom_prefix + chrom] = ft_inst.stop
            except KeyError:
                raise ValueError(
                    f"Feature {acc}: a region must have a 'Name' attribute.\n{print(ft_inst)}"
                )

        if get_genome_of_feature(ft_inst) in genome_skip or chrom in chrom_skip:
            continue
        elif not verify_feature_type(ft_inst):
            raise ValueError(f"Feature {acc}: unknown feature type '{ft_inst.featuretype}'.")
        elif ft_inst.featuretype in feature_types_to_skip:
            continue
        else:
            if chrom == "":
                raise AssertionError(f"No chromosome found for feature '{ft_inst.seqid}'.")

            # Skip gene features whose coordinates duplicate an RNA annotation
            if ft_inst.featuretype == 'gene' and ft_inst.id in genes_with_rna_children:
                continue

            ft_inst.id = get_feature_name(ft_inst, **kwargs)
            ft_inst.chrom = chrom_prefix + chrom
            try:
                specs = buffer_spec_index.get(ft_inst.id) or buffer_spec_index.get(ft_inst.featuretype, [])
                for spec in specs:
                    for bf in create_buffered_features(ft_inst,
                                                       buffer_spec=spec,
                                                       chrom_lengths=chrom_lengths):
                        bf.chrom = chrom_prefix + chrom
                        outf.write(format_feature_as_bed(bf) + "\n")
            except KeyError:
                pass
            except TypeError:
                raise TypeError(f"Feature '{ft_inst.id}' has an invalid type: {type(annot_features[ft_inst.id])}")
            outf.write(format_feature_as_bed(ft_inst) + "\n")
            # Replicate on opposite strand if needed
            if should_replicate_on_both_strands(ft_inst):
                opposite_feature = create_opposite_strand_feature(ft_inst)
                outf.write(format_feature_as_bed(opposite_feature) + "\n")

    if annot_features:
        for _ft_id, ft_list in annot_features.items():
            for ft_inst in ft_list:
                ft_inst.id = get_feature_name(ft_inst)
                ft_inst.chrom = chrom_id_map[ft_inst.chrom]
                if ft_inst.chrom in chrom_skip:
                    continue
                outf.write(format_feature_as_bed(ft_inst) + "\n")
                # Replicate on opposite strand if needed
                if should_replicate_on_both_strands(ft_inst):
                    opposite_feature = create_opposite_strand_feature(ft_inst)
                    outf.write(format_feature_as_bed(opposite_feature) + "\n")

    print(f"{acc} features written to {bed_file}")
    outf.close()
    return

def get_feature_name(feature: gffutils.Feature, **kwargs) -> str:
    if feature.featuretype == "gene":
        try:
            gene_name = feature.attributes.get("ID")[0]
            return kwargs.get("name_index", {})[gene_name]
        except KeyError:
            return f"{feature.id}(gene)"
    elif feature.featuretype == "sequence_feature":
        subtype = feature.attributes.get("gbkey")[0]
        notes = feature.attributes.get("Note")[0].split('; ')
        ft_name = notes[0]
        if '-' in ft_name:
            ft_name = ft_name.split('-')[0]
        try:
            return f"{ft_name}({subtype})"
        except ValueError:
            raise ValueError(f"Unable to parse sequence feature type: {str(feature)}")
    elif feature.featuretype == "mobile_genetic_element":
        ft_type = feature.attributes.get("mobile_element_type", [""])[0]
        subtype, _subtype_name = ft_type.split(":")
        for regex in MGE_RE_LIST:
            match = regex.match(ft_type)
            if match:
                return f"{match.group(2)}({subtype})"
        raise ValueError(f"Unable to parse mobile element type: {str(feature)}")
    elif feature.featuretype == "long_terminal_repeat":
        subtype = feature.attributes.get("Note")[0]
        regex = re.compile(r'(.*) LTR')
        if regex.match(subtype):
            return f"{regex.match(subtype).group(1)}(LTR)"
        else:
            return subtype
    elif feature.featuretype == "rRNA":
        subtype = feature.attributes.get("product")[0]
        regex = re.compile(r"(\d+S) ribosomal RNA")
        if regex.match(subtype):
            return f"{regex.match(subtype).group(1)} rRNA"
        else:
            return subtype
    elif feature.featuretype in FEATURE_TYPE_NAME_MAP:
        return FEATURE_TYPE_NAME_MAP[feature.featuretype]
    raise ValueError(f"Unable to retrieve name of feature type: {feature.featuretype}")

def remove_redundant_features(gff_db: gffutils.FeatureDB) -> None:
    prev = None
    merged_features = {}
    for feature in gff_db.all_features():
        
        if prev is None or feature.chrom != prev.chrom:
            prev = feature
            continue

        if feature.featuretype != prev.featuretype or feature.strand != prev.strand:
            key = (feature.chrom, feature.featuretype, feature.strand, feature.start, feature.stop)
        elif feature.start >= prev.start and feature.stop <= prev.stop:
            key = (prev.chrom, prev.featuretype, prev.strand, prev.start, prev.stop)
        else:
            key = (feature.chrom, feature.featuretype, feature.strand, feature.start, feature.stop)

        try:
            merged_features[key].append(feature)
            gff_db = gff_db.delete(feature)
        except KeyError:
            merged_features[key] = [feature]
            prev = feature

    # for key, features in merged_features.items():
    #     if len(features) == 1:
    #         continue
    #     gff_db = gff_db.update(features, merge_strategy="merge")
    return merged_features

def get_temporary_gff_file_path() -> str:
    return TMP_GFF_PATH

def write_features_to_file(features: List[gffutils.Feature], out_file: str = '') -> str:
    if out_file == '':
        out_file = get_temporary_gff_file_path()

    
    with open(out_file, "w") as outf:
        header = "##gff-version 3\n"
        outf.write(header)
        while features:
            ft = features.pop(0)
            outf.write(str(ft))
            if features:
                outf.write("\n")
            
    return out_file

def get_ribosomal_protein_genes(gff_db: gffutils.FeatureDB, subunit_index: Dict[int, str] = None) -> Dict[str, str]:
    gene_name_lookup = {}
    rps_regex = re.compile(r'ribosomal (\d+)S subunit protein (\w+)')
    if subunit_index is None:
        subunit_index = {
            40: 'S',
            60: 'L'
        }
    for ft in gff_db.all_features(featuretype='mRNA'):
        parent_gene = ft.attributes.get("Parent")[0]
        try:
            product = ft.attributes.get("product")[0]
        except TypeError:
            continue

        match = rps_regex.match(product)
        if match:
            try:
                subunit = subunit_index[int(match.group(1))]
            except KeyError:
                raise KeyError(f"Invalid ribosomal subunit: {match.group(1)}")
            gene_name_lookup[parent_gene] = f"RP{subunit}"
    return gene_name_lookup


def create_buffered_features(ft_inst: gffutils.Feature,
                             buffer_spec: BUFFER_SPEC = None,
                             **kwargs) -> List[gffutils.Feature]:
    """Add upstream and downstream buffer regions to the GFF database

    Args:
        gff_db (gffutils.FeatureDB): A GFF database loaded by gffutils.create_db
        buffer_features (list): A list of feature types to add buffer regions to
        buffer_length (int, optional): The length of the upstream and downstream features. Defaults to 1000.

    Raises:
        ValueError: If a feature has an invalid strand
    """
    new_features = []
    chrom_lengths = kwargs.get("chrom_lengths", {})

    # Load the lengths of the chromosomes if not provided
    if len(chrom_lengths) == 0:
        gff_db = kwargs.get("gff_db")
        for ft in gff_db.all_features(featuretype='region'):
            chrom_lengths[ft.chrom] = ft.stop

    if buffer_spec:
        prefixes = buffer_spec.prefix
        prefix_sep = buffer_spec.sep
        buffer_feature = buffer_spec.feature_name
        buffer_length = buffer_spec.length
    else:
        prefixes = kwargs.get('prefixes', ['upstream', 'downstream'])
        prefix_sep = kwargs.get('prefix_sep', '_')
        buffer_feature = kwargs.get("buffer_feature")
        buffer_length = kwargs.get("buffer_length")

    if type(prefixes) is str:
        raise AssertionError("Buffer prefix must be a list of length 1 or 2")

    if len(prefixes) == 1:
        prefixes = prefixes * 2
    elif len(prefixes) == 2:
        pass
    else:
        raise ValueError("Buffer prefix must be a list of length 1 or 2")

    for i in prefixes:
        VALID_FEATURE_TYPES.add(f"{i}{prefix_sep}{buffer_feature}")
        FEATURE_TYPE_NAME_MAP[f"{i}{prefix_sep}{buffer_feature}"] = f"{i}{prefix_sep}{buffer_feature}"

    if ft_inst.strand == "+":
        us_start = max(ft_inst.start - buffer_length, 1)
        us_stop = ft_inst.start # - 1
        ds_start = ft_inst.stop # + 1
        ds_stop = min(ft_inst.stop + buffer_length, chrom_lengths[ft_inst.chrom])
    elif ft_inst.strand == "-":
        us_start = ft_inst.stop
        us_stop = min(ft_inst.stop + buffer_length, chrom_lengths[ft_inst.chrom])
        ds_start = max(ft_inst.start - buffer_length, 1)
        ds_stop = ft_inst.start
    else:
        raise ValueError(f"Feature {ft_inst.id} has an invalid strand: {ft_inst.strand}")

    for prefix, coords in zip(prefixes, [[us_start, us_stop], [ds_start, ds_stop]]):
        ft_type = f"{prefix}{prefix_sep}{ft_inst.featuretype}"
        ft_id = f"{ft_inst.id}{prefix_sep}{prefix}"
        buff_feat = feature_mod.feature_from_line(
                f"{ft_inst.chrom}\t{ft_inst.source}\t{ft_type}\t{coords[0]}\t{coords[1]}\t.\t{ft_inst.strand}\t.\tID={ft_id}"
            )
        buff_feat.id = ft_id
        new_features.append(buff_feat)
    
    # Add the new features to the GFF database
    # NOTE: COULD NOT GET THIS TO WORK
    # Option 1: Write features to file, then load database and update
    # all_new_features = [item for sublist in new_features.values() for item in sublist]
    # new_features_gff = write_features_to_file(all_new_features)
    # new_gff = create_gff_db(new_features_gff)
    # gff_db = gff_db.update(new_gff)

    # Option 2: Update the database directly
    # gff_db = gff_db.update(new_features)

    return new_features

def load_buffer_feature_specifications(buffer_spec_file: str) -> List[BUFFER_SPEC]:
    if not os.path.exists(buffer_spec_file):
        raise FileNotFoundError(f"Buffer specification file not found: {buffer_spec_file}")
    
    buffer_specs = []

    with open(buffer_spec_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                try:
                    ft_name, prefix, prefix_sep, length = line.strip().split("\t")
                except ValueError:
                    raise ValueError(f"Invalid buffer specification: {line}")

                try:
                    buffer_specs.append(BUFFER_SPEC(
                        feature_name=ft_name,
                        prefix=prefix.split(','),
                        sep=prefix_sep,
                        length=int(length)
                        ))
                except ValueError:
                    raise ValueError(f"Unable to use feature buffer specification:\n{line}\n")

    return buffer_specs

def load_annotations_from_bed(annot_table: str) -> Dict[str, list]:
    if not os.path.exists(annot_table):
        raise FileNotFoundError(f"Annotation table not found: {annot_table}")

    annot_dict = {}
    i = 0
    with open(annot_table) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            chrom, start, end, ft_name, strand = line.strip().split("\t")
            annot_dict[i] = [chrom, start, end, ft_name, strand]
            i += 1

    return annot_dict


def get_buffer_specification_index(buffer_spec_file: str = '',
                                   **kwargs) -> Dict[str, List[BUFFER_SPEC]]:
    buffer_spec_dict = {}
    buffer_specs = kwargs.get("buffer_spec", {})
    
    if not buffer_specs:  
        if not buffer_spec_file:
            return

        buffer_specs = load_buffer_feature_specifications(buffer_spec_file)
    
    for spec in buffer_specs:
            try:
                buffer_spec_dict[spec.feature_name].append(spec)
            except KeyError:
                buffer_spec_dict[spec.feature_name] = [spec]
    return buffer_spec_dict

def create_subtelomere_features(gff_db: gffutils.FeatureDB, **kwargs) -> Dict[str, List[gffutils.Feature]]:
    subtelomere_features = {}
    feature_length = kwargs.get('length', 15000)
    for ft_inst in gff_db.all_features(featuretype='telomere'):
        # Create a unique name for the subtelomere feature
        subtelomere_name = f"subtelomere-{ft_inst.chrom}"
        if subtelomere_name not in subtelomere_features:
            subtelomere_features[subtelomere_name] = []

        if ft_inst.start == 1 or ft_inst.stop == 1:
            # This is the 5' subtelomere
            start = max(ft_inst.start, ft_inst.stop)
            stop = start + feature_length
        else:
            # This is the 3' subtelomere
            stop = min(ft_inst.start, ft_inst.stop)
            start = stop - feature_length

        if kwargs.get('flip_negative_strand_positions', False):
            if ft_inst.strand == "+":
                pass
            elif ft_inst.strand == "-":
                start, stop = stop, start
            else:
                raise ValueError(f"Invalid strand: {ft_inst.strand}")

        subtelomere_features[subtelomere_name] += [
                feature_mod.feature_from_line(
                    f"{ft_inst.chrom}\t{ft_inst.source}\tsubtelomere\t{start}\t{stop}\t.\t{ft_inst.strand}\t.\tID={subtelomere_name}-{len(subtelomere_features[subtelomere_name])}"
                )
            ]
        for subt_inst in subtelomere_features[subtelomere_name]:
            subt_inst.id = subtelomere_name
            subt_inst.featuretype = "subtelomere"
    return subtelomere_features

def update_gff_with_additional_annotations(gff_db: gffutils.FeatureDB,
                                           annot_table: str = '') ->  Dict[str, List[gffutils.Feature]]:
    """Add user-specified annotations to the GFF database

    This function can be used to add non-canonical genome feature annotations
    such as the pericentromere, subtelomere, or other regions.

    The annot_table should be a BED file with the following columns:
    1. Chromosome
    2. Start position
    3. End position
    4. Feature name
    5. Strand

    The chromosome can use the special values 'all' or 'mitochondrion' to apply
    the annotation to all nuclear chromosomes or the mitochondrion, respectively.
    Otherwise, the chromosome should match the chromosome names in the GFF file.

    Args:
        gff_db (gffutils.FeatureDB): A GFF database loaded by gffutils.create_db
        annot_table (str): Path to a BED file containing additional annotations
         that were not included in the original GFF file
    """
    new_features = {}

    def _get_chromosome_names(gff_db: gffutils.FeatureDB) -> List[str]:
        return [ft.chrom for ft in gff_db.all_features(featuretype='region')]
    
    def _create_feature_from_line(name: str, chrom: str, start: int, stop: int, strand: str, **kwargs) -> gffutils.Feature:
        ft_type = kwargs.get('feature_type', '')
        if ft_type and ft_type not in VALID_FEATURE_TYPES:
            raise ValueError(f"Invalid feature type: {ft_type}")
        elif ft_type != '':
            ft_id = f"{ft_type}-{name}"
        else:
            ft_id = name

        ft = feature_mod.feature_from_line(
            f"{chrom}\t.\t{name}\t{start}\t{stop}\t.\t{strand}\t.\tID={ft_id}"
        )
        ft.id = ft_id
        ft.featuretype = ft_type
        return ft
    
    if annot_table:
        annot_dict = load_annotations_from_bed(annot_table)
    else:
        annot_dict = {}
    
    chrom_names = _get_chromosome_names(gff_db)
    for i, ft_coords in annot_dict.items():
        chrom, start, stop, ft_name, strand = ft_coords
        if chrom not in chrom_names and chrom not in ['all', 'mitochondrion']:
            raise ValueError(f"Chromosome '{chrom}' not found in GFF file")

        if chrom == 'all':
            for chrom in chrom_names:
                ft = _create_feature_from_line(
                    f"{ft_name}-{i}", chrom, int(start), int(stop), strand, feature_type=ft_name
                )
                # gff_db.add_feature(ft)
                try:
                    new_features[ft.id].append(ft)
                except KeyError:
                    new_features[ft.id] = [ft]
        else:
            ft = _create_feature_from_line(ft_name, chrom, int(start), int(stop), strand)
            # gff_db.add_feature(ft)
            new_features[ft.id] = [ft]

    return new_features

def main(args=None):
    opts = get_options(args)
    print("Loading GFF into database")
    gff_db = create_gff_db(opts.p_i_gff)
    print("Merging redundant features")
    remove_redundant_features(gff_db)
    print("Adding buffer regions")
    buffer_features = get_buffer_specification_index(opts.p_i_buffer_features)
    
    print("Adding additional annotations")
    annot_features = update_gff_with_additional_annotations(gff_db,
                                                            annot_table=opts.p_i_annot)
    annot_features.update(create_subtelomere_features(gff_db))

    gene_name_lookup = get_ribosomal_protein_genes(gff_db)

    print("Writing GFF features to BED")
    parse_gff_to_bed(gff_db,
                     opts.p_o_bed,
                     annot_features=annot_features,
                     buffer_spec_index=buffer_features,
                     name_index=gene_name_lookup)
    return


if __name__ == "__main__":
    main()
