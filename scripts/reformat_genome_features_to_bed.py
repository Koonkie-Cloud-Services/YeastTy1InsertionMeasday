#!/usr/bin/env python

import argparse
import sys
from typing import List, Dict

import gffutils
from gffutils import inspect
from gffutils import feature
from gffutils import create

FEATURE_TYPE_NAME_MAP = {
    "gene": "gene",
    "telomere": "telomere",
    "centromere": "centromere",
    "origin_of_replication": "rep_origin",
    "tRNA": "tRNA",
    "ncRNA": "ncRNA",
    "snoRNA": "snoRNA",
    "snRNA": "snRNA",
    "rRNA": "rRNA",
    "SRP_RNA": "SRP_RNA",
    "pseudogene": "pseudogene",
    "antisense_RNA": "antisense_RNA",
    "telomerase_RNA": "telomerase_RNA",
    "sequence_feature": "sequence_feature",
    "mobile_genetic_element": "mobile_element",
    "repeat_region": "LTR",
    "long_terminal_repeat": "LTR",
}
SKIPPED_FEATURE_TYPES = [
    "mRNA",
    "CDS",
    "exon",
    "intron",
    "region",
    "transcript",
    "RNase_P_RNA",
    "RNase_MRP_RNA",
    "regulatory_region",
]
VALID_FEATURE_TYPES = list(FEATURE_TYPE_NAME_MAP) + SKIPPED_FEATURE_TYPES
TMP_GFF_PATH = "tmp.gff"

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
        dest="p_buffer_features",
        nargs="+",
        default=["tRNA"],
        help="A list of features to add surrounding upstream and downstream region features"
    )
    parser.add_argument(
        "--buffer_length",
        dest="p_buffer_length",
        type=int,
        default=1000,
        help="The length of the buffer regions to add"
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
        get_feature_name(feature),
        feature.strand
    ]
    return "\t".join([str(i) for i in bed_row])

def parse_gff_to_bed(gff_db: gffutils.FeatureDB, bed_file: str, **kwargs) -> None:
    outf = open(bed_file, "w")

    chrom_prefix = kwargs.get("chrom_prefix", "chr")
    genome_skip = kwargs.get("genome_skip", ['mitochondrion', 'plasmid'])
    chrom_skip = kwargs.get("chrom_skip", ['MT'])
    feature_types_to_skip = kwargs.get("feature_types_to_skip", SKIPPED_FEATURE_TYPES)

    buffer_regions = kwargs.get("buffer_regions", {})
    annot_features = kwargs.get("annot_features", {})

    chrom = ""
    acc = 0
    for ft_inst in gff_db.all_features():
        acc += 1

        if ft_inst.featuretype == "region":
            # Set the chromosome name for the following features
            try:
                chrom = ft_inst.attributes.get("Name")[0]
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

            ft_inst.chrom = chrom_prefix + chrom
            try:
                buffer_ft = buffer_regions[ft_inst.id]
                for ft in buffer_ft:
                    ft.chrom = chrom_prefix + chrom
                    outf.write(format_feature_as_bed(ft) + "\n")
            except KeyError:
                pass
            outf.write(format_feature_as_bed(ft_inst) + "\n")

    print(f"{acc} features written to {bed_file}")
    outf.close()
    return

def get_feature_name(feature: gffutils.Feature) -> str:
    if feature.featuretype == "gene":
        return f"{feature.id}(gene)"
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

def add_buffer_regions(gff_db: gffutils.FeatureDB, buffer_features: list, buffer_length = 1000) -> Dict[str, List[gffutils.Feature]]:
    """Add upstream and downstream buffer regions to the GFF database

    Args:
        gff_db (gffutils.FeatureDB): A GFF database loaded by gffutils.create_db
        buffer_features (list): A list of feature types to add buffer regions to
        buffer_length (int, optional): The length of the upstream and downstream features. Defaults to 1000.

    Raises:
        ValueError: If a feature has an invalid strand
    """
    new_features = {}
    chrom_lengths = {}

    for ft_name in buffer_features:
        VALID_FEATURE_TYPES.append(f"upstream_{ft_name}")
        VALID_FEATURE_TYPES.append(f"downstream_{ft_name}")
        FEATURE_TYPE_NAME_MAP[f"upstream_{ft_name}"] = f"upstream_{ft_name}"
        FEATURE_TYPE_NAME_MAP[f"downstream_{ft_name}"] = f"downstream_{ft_name}"

    for ft_inst in gff_db.all_features():
        if ft_inst.featuretype == "region":
            # get chromosome length
            chrom_lengths[ft_inst.chrom] = ft_inst.stop
        if ft_inst.featuretype in buffer_features:
            if ft_inst.strand == "+":
                us_start = max(ft_inst.start - buffer_length, 0)
                us_stop = ft_inst.start # - 1
                ds_start = ft_inst.stop # + 1
                ds_stop = min(ft_inst.stop + buffer_length, chrom_lengths[ft_inst.chrom])
            elif ft_inst.strand == "-":
                us_start = min(ft_inst.stop, chrom_lengths[ft_inst.chrom])
                us_stop = ft_inst.stop + buffer_length # + 1
                ds_start = ft_inst.start - buffer_length # - 1
                ds_stop = max(ft_inst.start, 0)
            else:
                raise ValueError(f"Feature {ft_inst.id} has an invalid strand: {ft_inst.strand}")
        
            upstream = feature.feature_from_line(
                        f"{ft_inst.chrom}\t{ft_inst.source}\tupstream_{ft_inst.featuretype}\t{us_start}\t{us_stop}\t.\t{ft_inst.strand}\t.\tID={ft_inst.id}_upstream"
                    )
            upstream.id = f"{ft_inst.id}_upstream"
            downstream = feature.feature_from_line(
                f"{ft_inst.chrom}\t{ft_inst.source}\tdownstream_{ft_inst.featuretype}\t{ds_start}\t{ds_stop}\t.\t{ft_inst.strand}\t.\tID={ft_inst.id}_downstream"
            )
            downstream.id = f"{ft_inst.id}_downstream"
            new_features[ft_inst.id] = [upstream, downstream]
        else:
            pass
    
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

def update_gff_with_additional_annotations(gff_db: gffutils.FeatureDB, annot_table: str) -> None:
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
    return

def main():
    opts = get_options()
    print("Loading GFF into database")
    gff_db = create_gff_db(opts.p_i_gff)
    print("Merging redundant features")
    remove_redundant_features(gff_db)
    print("Adding buffer regions")
    buffer_regions = add_buffer_regions(gff_db,
                                        buffer_features=opts.p_buffer_features,
                                        buffer_length=opts.p_buffer_length)
    print("Adding additional annotations")
    annot_features = update_gff_with_additional_annotations(gff_db,
                                                            annot_table=opts.p_i_annot)
    print("Parsing GFF to BED")
    parse_gff_to_bed(gff_db,
                     opts.p_o_bed,
                     buffer_regions=buffer_regions,
                     annot_features=annot_features)
    return


if __name__ == "__main__":
    main()
