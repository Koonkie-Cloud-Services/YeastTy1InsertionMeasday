#!/usr/bin/env python

import gffutils
import argparse
import sys


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

def parse_gff_to_bed(gff_db: gffutils.FeatureDB, bed_file: str, **kwargs) -> None:
    outf = open(bed_file, "w")

    chrom_prefix = kwargs.get("chrom_prefix", "chr")
    genome_skip = kwargs.get("genome_skip", ['mitochondrion', 'plasmid'])
    chrom = ""
    acc = 0
    for feature in gff_db.all_features():
        acc += 1
        if get_genome_of_feature(feature) in genome_skip:
            continue

        if feature.featuretype in [
            "gene",
            "telomere",
            "centromere",
            "origin_of_replication",
            "tRNA",
            "repeat_region",
            "mobile_element",
        ]:
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
        elif feature.featuretype == "region":
            try:
                chrom = feature.attributes["chromosome"][0]
            except KeyError:
                raise ValueError(
                    f"Feature {acc}: a region must have a chromosome attribute.\n{print(feature)}"
                )
        else:
            continue

    outf.close()
    return


def get_feature_name(feature: gffutils.Feature) -> str:
    feature_type_map = {
        "telomere": "telomere",
        "centromere": "centromere",
        "origin_of_replication": "rep_origin",
        "tRNA": "tRNA",
        "repeat_region": "LTR",
        "mobile_element": "LTR",
    }
    if feature.featuretype == "gene":
        return f"{feature.id}(gene)"
    elif feature.featuretype in feature_type_map:
        return feature_type_map[feature.featuretype]
    raise ValueError(f"Unknown feature type: {feature.featuretype}")


def main():
    opts = get_options()
    gff_db = create_gff_db(opts.p_i_gff)
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
