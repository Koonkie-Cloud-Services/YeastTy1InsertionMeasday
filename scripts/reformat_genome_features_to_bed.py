#!/usr/bin/env python

import gffutils
import argparse
import sys

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


def create_gff_db(gff_file) -> gffutils.FeatureDB:
    gff_db = gffutils.create_db(
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
    chrom = ""
    acc = 0
    for feature in gff_db.all_features():
        acc += 1

        if feature.featuretype == "region":
            # Set the chromosome name for the following features
            try:
                chrom = feature.attributes.get("Name")[0]
            except KeyError:
                raise ValueError(
                    f"Feature {acc}: a region must have a 'Name' attribute.\n{print(feature)}"
                )

        if get_genome_of_feature(feature) in genome_skip or chrom in chrom_skip:
            continue
        elif not verify_feature_type(feature):
            raise ValueError(f"Feature {acc}: unknown feature type '{feature.featuretype}'.")
        elif feature.featuretype in feature_types_to_skip:
            continue
        else:
            if chrom == "":
                raise AssertionError(f"No chromosome found for feature '{feature.seqid}'.")

            feature.chrom = chrom_prefix + chrom
            outf.write(format_feature_as_bed(feature) + "\n")

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

def add_buffer_regions(gff_db: gffutils.FeatureDB, buffer_features: list, buffer_length = 1000) -> None:
    return

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
    add_buffer_regions(gff_db,
                       buffer_features=opts.p_buffer_features,
                       buffer_length=opts.p_buffer_length)
    print("Adding additional annotations")
    update_gff_with_additional_annotations(gff_db,
                                           annot_table=opts.p_i_annot)
    print("Parsing GFF to BED")
    parse_gff_to_bed(gff_db, opts.p_o_bed)
    return


if __name__ == "__main__":
    main()
