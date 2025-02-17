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
        "--bed_out",
        dest="p_o_bed",
        type=str,
        help="The name of the BED file to create"
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
            'featuretype',
            'score',
            'strand',
            'frame'
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

def parse_gff_to_bed(gff_db: gffutils.FeatureDB, bed_file: str, **kwargs) -> None:
    outf = open(bed_file, "w")

    chrom_prefix = kwargs.get("chrom_prefix", "chr")
    genome_skip = kwargs.get("genome_skip", ['mitochondrion', 'plasmid'])
    chrom_skip = kwargs.get("chrom_skip", ['MT'])
    chrom = ""
    acc = 0
    for feature in gff_db.all_features():
        acc += 1

        if feature.featuretype == "region":
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
        elif feature.featuretype in SKIPPED_FEATURE_TYPES:
            continue
        else:
            if chrom == "":
                raise AssertionError(f"No chromosome found for feature '{feature.seqid}'.")

            bed_row = [
                chrom_prefix + chrom,
                feature.start,
                feature.stop,
                get_feature_name(feature),
                feature.strand
            ]
            outf.write(
                "\t".join([str(i) for i in bed_row]) + "\n"
            )

    outf.close()
    return


def get_feature_name(feature: gffutils.Feature) -> str:
    if feature.featuretype == "gene":
        return f"{feature.id}(gene)"
    elif feature.featuretype in FEATURE_TYPE_NAME_MAP:
        return FEATURE_TYPE_NAME_MAP[feature.featuretype]
    raise ValueError(f"Unable to retrieve name of feature type: {feature.featuretype}")


def main():
    opts = get_options()
    gff_db = create_gff_db(opts.p_i_gff)
    parse_gff_to_bed(gff_db, opts.p_o_bed)
    return


if __name__ == "__main__":
    main()
