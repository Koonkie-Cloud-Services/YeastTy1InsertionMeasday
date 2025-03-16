#!/usr/bin/env python

import argparse
import sys
from typing import List, Dict
import pandas as pd
from gffutils import FeatureDB as gff_FeatureDB
from gffutils import Feature as gff_Feature

from reformat_genome_features_to_bed import create_gff_db
from reformat_genome_features_to_bed import BUFFER_SPEC


class Feature(gff_Feature):
    def __init__(
        self,
        featuretype: str,
        seqid: str = ".",
        start: int = -1,
        end: int = -1,
        strand: str = "",
        **kwargs
    ):
        self.seqid = str(seqid)
        self.start = int(start)
        self.end = int(end)
        self.featuretype = str(featuretype)
        self.strand = str(strand)
        self.attributes = {}

        if len(kwargs) > 0:
            self.attributes[kwargs.get("Attribute_name")] = kwargs.get("Attribute_value")
    
        return

    def __str__(self):
        return f"{self.seqid}\t{self.start}\t{self.end}\t{self.featuretype}\t{self.strand}"

    def has_coordinates(self) -> bool:
        return self.start >= 0 and self.end > 1 and self.seqid != ""


class FeatureDB:
    def __init__(self, features: List[Feature]):
        self.features = features

    def __iter__(self):
        self._index = 0
        return self

    def __next__(self):
        if self._index < len(self.features):
            result = self.features[self._index].featuretype
            self._index += 1
            return result
        else:
            raise StopIteration

    def get_features_of_type(self, feature_type: str) -> List[Feature]:
        return [feature for feature in self.features if feature.featuretype == feature_type]
    
    def all_features(self, **kwargs) -> List[Feature]:
        if kwargs.get("featuretype") is None:
            return self.features
        return self.get_features_of_type(kwargs.get("featuretype"))


def feature_from_bed_line(line: str) -> Feature:
    seqid, start, end, featuretype, strand = line.strip().split("\t")
    return Feature(featuretype, seqid=seqid, start=start, end=end, strand=strand)


def get_options(args=sys.argv[1:]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert the features in a GFF file to BED format"
    )
    features_input_parser = parser.add_mutually_exclusive_group(required=True)
    features_input_parser.add_argument(
        "--bed_in",
        dest="p_i_bed",
        type=str,
        help="Path to the input BED file containing genomic features",
    )
    features_input_parser.add_argument(
        "--gff_in",
        dest="p_i_gff",
        type=str,
        help="Path to the input GFF file containing genomic features"
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
        default=0,
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
        if not pd.isna(row["Attribute.name"]):
            feature = Feature(
                featuretype=row["Feature.name"],
                Attribute_name=row["Attribute.name"],
                Attribute_value=row["Attribute.value"],
            )
        elif pd.isna(row["ID.name"]):
            feature = Feature(
                featuretype=row["Feature.name"],
            )
        else:
            feature = Feature(
                seqid=row["Chromosome"],
                strand=row["Strand"],
                start=int(row["Start"]),
                end=int(row["Stop"]),
                featuretype=row["Feature.name"],
            )
        try:
            target_features[row["RNA.Polymerase"]].append(feature)
        except KeyError:
            target_features[row["RNA.Polymerase"]] = [feature]
    return target_features


def match_features_in_db(
    feature_db: FeatureDB | gff_FeatureDB, features_in: List[Feature]
) -> List[Feature]:
    features_out = []
    for feature in features_in:
        if feature.has_coordinates():
            features_out.append(feature)
        else:
            for db_feat in feature_db.all_features(featuretype=feature.featuretype):
                if len(feature.attributes):
                    for attr_name, attr_value in feature.attributes.items():
                        if db_feat.attributes.get(attr_name)[0] == attr_value:
                            features_out.append(db_feat)
                else:
                    features_out.append(db_feat)
    return features_out


def load_buffer_specs(p_i_xlsx: str,
                      p_i_sheet_name: str,
                      **kwargs) -> Dict[str, BUFFER_SPEC]:
    buffer_specs = {}
    df = pd.read_excel(p_i_xlsx, sheet_name=p_i_sheet_name)
    for _index, row in df.iterrows():
        if pd.isna(row["Buffer.length"]):
            bl = kwargs.get("default_buffer_length_bp", 0)
        else:
            bl = row["Buffer.length"]

        buffer_specs[row["Feature.name"]] = BUFFER_SPEC(
            feature_name = row["Feature.name"],
            prefix = "",
            sep = " ",
            length = bl
        )
    return buffer_specs


def buffer_features(features_in: List[Feature],
                    buffer_specs: Dict[str, BUFFER_SPEC]) -> List[Feature]:
    features_out = []
    for feature in features_in:
        buffer_spec = buffer_specs.get(feature.featuretype)
        if buffer_spec is None:
            features_out.append(feature)
            continue

        us_feat = Feature(
            featuretype=f"upstream {feature.featuretype}",
            seqid=feature.seqid,
            start=max(feature.start - buffer_spec.length, 1),
            end=feature.start,
            strand=feature.strand,
        )
        ds_feat = Feature(
            featuretype=f"downstream {feature.featuretype}",
            seqid=feature.seqid,
            start=feature.end,
            end=feature.end + buffer_spec.length,
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
    if opts.p_i_gff:
        feature_db = create_gff_db(opts.p_i_gff)
    elif opts.p_i_bed:
        feature_db = load_features_from_bed_file(opts.p_i_bed)
    else:
        raise ValueError("No input file provided")

    target_features = load_target_feature_info(opts.p_i_xlsx,
                                               opts.p_i_sheet_name)
    buffer_specs = load_buffer_specs(opts.p_i_xlsx,
                                     opts.p_i_sheet_name,
                                     default_buffer_length_bp=1000)

    for rna_pol, feature_set in target_features.items():
        output_bed_path = get_bed_file_name(opts.p_o_bed_prefix, rna_pol)
        feature_set = match_features_in_db(feature_db, feature_set)
        feature_set = buffer_features(feature_set, buffer_specs)
        write_features_to_bed_file(output_bed_path, feature_set)
    return


if __name__ == "__main__":
    main()
