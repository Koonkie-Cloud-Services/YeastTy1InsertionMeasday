#!/usr/bin/env python

import argparse
import sys
from typing import List, Dict
import pandas as pd


class Feature:
    def __init__(
        self,
        name: str,
        chrom: str = "",
        start: int = -1,
        end: int = -1,
        strand: str = "",
    ):
        self.chrom = str(chrom)
        self.start = int(start)
        self.end = int(end)
        self.name = str(name)
        self.strand = str(strand)

    def __str__(self):
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.name}\t{self.strand}"

    def has_coordinates(self) -> bool:
        return self.start >= 0 and self.end > 1 and self.chrom != ""


class FeatureDB:
    def __init__(self, features: List[Feature]):
        self.features = features

    def __iter__(self):
        self._index = 0
        return self

    def __next__(self):
        if self._index < len(self.features):
            result = self.features[self._index].name
            self._index += 1
            return result
        else:
            raise StopIteration

    def get_features_of_type(self, feature_type: str) -> List[Feature]:
        return [feature for feature in self.features if feature.name == feature_type]


def feature_from_bed_line(line: str) -> Feature:
    chrom, start, end, name, strand = line.strip().split("\t")
    return Feature(name, chrom=chrom, start=start, end=end, strand=strand)


def get_options(args=sys.argv[1:]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert the features in a GFF file to BED format"
    )
    parser.add_argument(
        "--bed_in",
        dest="p_i_bed",
        type=str,
        help="Path to the input BED file containing genomic features",
    )
    parser.add_argument(
        "--xlsx_in",
        dest="p_i_xlsx",
        type=str,
        required=False,
        help="An Excel file (.xlsx) containing the features to write to a BED file.",
    )
    parser.add_argument(
        "--bed_out_prefix",
        dest="p_o_bed_prefix",
        type=str,
        help="The prefix of the BED files to create."
        "RNA polymerase will be appended to the prefix.",
    )
    parser.add_argument(
        "--buffer_length_bp",
        dest="p_i_buffer_features",
        type=int,
        default=1000,
        required=False,
        help="The number of base pairs to add to the start and end of each feature",
    )
    parser.add_argument(
        "--sheet_name",
        dest="p_i_sheet_name",
        type=str,
        required=True,
        help="The name of the sheet in xlsx_in containing feature information",
    )
    return parser.parse_args(args=args)


def get_bed_file_name(p_o_bed_prefix: str, rna_pol: int) -> str:
    return f"{p_o_bed_prefix}.RNA_Pol_{rna_pol}.bed"

def load_features_from_bed_file(p_i_bed: str) -> FeatureDB:
    features = []
    with open(p_i_bed, "r") as f:
        for line in f:
            features.append(feature_from_bed_line(line))
    return FeatureDB(features)


def load_target_feature_info(
    p_i_xlsx: str, p_i_sheet_name: str
) -> Dict[str, List[Feature]]:
    target_features = {}
    df = pd.read_excel(p_i_xlsx, sheet_name=p_i_sheet_name)
    for index, row in df.iterrows():
        if pd.isna(row["ID.name"]):
            feature = Feature(
                name=row["Feature.name"],
            )
        else:
            feature = Feature(
                chrom=row["Chromosome"],
                strand=row["Strand"],
                start=int(row["Start"]),
                end=int(row["Stop"]),
                name=row["Feature.name"],
            )
        try:
            target_features[row["RNA.Polymerase"]].append(feature)
        except KeyError:
            target_features[row["RNA.Polymerase"]] = [feature]
    return target_features


def match_features_in_db(
    feature_db: FeatureDB, features_in: List[Feature]
) -> List[Feature]:
    features_out = []
    for feature in features_in:
        if feature.has_coordinates():
            features_out.append(feature)
        else:
            for db_feat in feature_db.get_features_of_type(feature.name):
                features_out.append(db_feat)
    return features_out


def buffer_features(features_in: List[Feature]) -> List[Feature]:
    features_out = []
    for feature in features_in:
        us_feat = Feature(
            name=f"upstream {feature.name}",
            chrom=feature.chrom,
            start=max(feature.start - 1000, 1),
            end=feature.start,
            strand=feature.strand,
        )
        ds_feat = Feature(
            name=f"downstream {feature.name}",
            chrom=feature.chrom,
            start=feature.end,
            end=feature.end + 1000,
            strand=feature.strand,
        )
        features_out.extend([us_feat, feature, ds_feat])
    return features_out


def write_features_to_bed_file(p_o_bed: str, features_in: List[Feature]) -> None:
    with open(p_o_bed, "w") as f:
        for feature in features_in:
            f.write(f"{str(feature)}\n")
    return


def main(args=None):
    opts = get_options(args)
    feature_db = load_features_from_bed_file(opts.p_i_bed)
    target_features = load_target_feature_info(opts.p_i_xlsx, opts.p_i_sheet_name)
    for rna_pol, feature_set in target_features.items():
        output_bed_path = get_bed_file_name(opts.p_o_bed_prefix, rna_pol)
        feature_set = match_features_in_db(feature_db, feature_set)
        feature_set = buffer_features(feature_set)
        write_features_to_bed_file(output_bed_path, feature_set)
    return


if __name__ == "__main__":
    main()
