#!/usr/bin/env python

import gffutils
import argparse
import sys


def get_options(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        description="Convert the features in a GFF file to BED format"
    )
    parser.add_argument(
        "--p_i_gff", type="str", help="The GFF file to convert"
    )
    parser.add_argument(
        "--p_o_bed", type="str", help="The name of the BED file to create"
    )
    return parser.parse_args(args=args)


def create_gff_db(gff_file) -> str:
    db_fn = gff_file + ".db"
    gffutils.create_db(
        gff_file,
        db_fn,
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
    )
    return db_fn


def parse_gff_to_bed(gff_db: gffutils.FeatureDB, bed_file: str) -> None: 
    outf = open(bed_file, "w")
    for feature in gff_db.all_features():
        if feature.featuretype in [
            "gene",
            "telomere",
            "centromere",
            "rep_origin",
            "tRNA",
            "repeat_region",
            "mobile_element",
        ]:
            name = get_feature_name(feature)
            bed_line = "{0}\t{1}\t{2}\t{3}\n".format(feature.start, feature.stop, name, feature.strand)
            outf.write(bed_line)
    outf.close()
    return


def get_feature_name(feature: gffutils.Feature) -> str:
    if feature.featuretype == "gene":
        return f"{feature.id}(gene)"
    elif feature.featuretype in ["telomere", "centromere", "rep_origin", "tRNA"]:
        return feature.featuretype
    elif feature.featuretype in ["repeat_region", "mobile_element"]:
        return "LTR"
    raise ValueError(f"Unknown feature type: {feature.featuretype}")


def main():
    opts = get_options()
    gff_db = gffutils.FeatureDB(create_gff_db(opts.p_i_gff))
    parse_gff_to_bed(gff_db, opts.p_o_bed)
    return


    # for record in SeqIO.parse(open(sys.argv[1], "rU"), "genbank"):
    #     for feature in record.features:
    #         if feature.type == "gene":
    #             start = feature.location.start.position
    #             stop = feature.location.end.position
    #             try:
    #                 # name = "gene"
    #                 name = feature.qualifiers["gene"][0] + "(gene)"
    #             except:
    #                 # some features only have a locus tag
    #                 name = feature.qualifiers["locus_tag"][0] + "(gene)"
    #             if feature.strand < 0:
    #                 strand = "-"
    #             else:
    #                 strand = "+"
    #             bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
    #             outf.write(bed_line)
    #         if feature.type == "telomere":
    #             start = feature.location.start.position
    #             stop = feature.location.end.position
    #             try:
    #                 name = "telomere"
    #             except:
    #                 name = "false"
    #             if feature.strand < 0:
    #                 strand = "-"
    #             else:
    #                 strand = "+"
    #             bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
    #             outf.write(bed_line)
    #         if feature.type == "centromere":
    #             start = feature.location.start.position
    #             stop = feature.location.end.position
    #             try:
    #                 name = "centromere"
    #             except:
    #                 name = "false"
    #             if feature.strand < 0:
    #                 strand = "-"
    #             else:
    #                 strand = "+"
    #             bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
    #             outf.write(bed_line)
    #         if feature.type == "rep_origin":
    #             start = feature.location.start.position
    #             stop = feature.location.end.position
    #             try:
    #                 name = "rep_origin"
    #             except:
    #                 name = "false"
    #             if feature.strand < 0:
    #                 strand = "-"
    #             else:
    #                 strand = "+"
    #             bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
    #             outf.write(bed_line)
    #         if feature.type == "tRNA":
    #             start = feature.location.start.position
    #             stop = feature.location.end.position
    #             try:
    #                 name = "tRNA"
    #             except:
    #                 name = "false"
    #             if feature.strand < 0:
    #                 strand = "-"
    #             else:
    #                 strand = "+"
    #             bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
    #             outf.write(bed_line)
    #         if feature.type == "repeat_region":
    #             start = feature.location.start.position
    #             stop = feature.location.end.position
    #             try:
    #                 name = "LTR"
    #             except:
    #                 name = "false"
    #             if feature.strand < 0:
    #                 strand = "-"
    #             else:
    #                 strand = "+"
    #             bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
    #             outf.write(bed_line)
    #         if feature.type == "mobile_element":
    #             start = feature.location.start.position
    #             stop = feature.location.end.position
    #             try:
    #                 name = "LTR"
    #             except:
    #                 name = "false"
    #             if feature.strand < 0:
    #                 strand = "-"
    #             else:
    #                 strand = "+"
    #             bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
    #             outf.write(bed_line)
    #     outf.close()


if __name__ == "__main__":
    main()
