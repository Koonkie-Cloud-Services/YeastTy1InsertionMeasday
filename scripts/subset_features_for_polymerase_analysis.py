#!/usr/bin/env python

import argparse
import sys
import re
from typing import List, Dict
import pandas as pd
from gffutils import FeatureDB as gff_FeatureDB
from gffutils import Feature as gff_Feature

from reformat_genome_features_to_bed import create_gff_db
from reformat_genome_features_to_bed import BUFFER_SPEC
from reformat_genome_features_to_bed import format_feature_as_bed


# Define global settings for each RNA polymerase
RNA_POL_SETTINGS = {
    3: {
        "feature_attribute": "id",
    },
    2: {
        "feature_attribute": "featuretype",
    },
}

class Feature(gff_Feature):
    def __init__(
        self,
        featuretype: str,
        id: str = "",
        seqid: str = ".",
        start: int = -1,
        end: int = -1,
        strand: str = "",
        **kwargs,
    ):
        self.seqid = str(seqid)
        self.start = int(start)
        self.end = int(end)
        self.featuretype = str(featuretype)
        self.strand = str(strand)
        self.id = id
        self.attributes = {}

        if len(kwargs) > 0:
            self.attributes[kwargs.get("Attribute_name")] = kwargs.get(
                "Attribute_value"
            )

        return

    def to_bed(self, attr_name="featuretype") -> str:
        try:
            ft_name = self.__getattribute__(attr_name)
        except AttributeError:
            raise AttributeError(
                f"Feature {self} does not have the attribute {attr_name}"
            )

        # Make sure the attribute value is valid
        if ft_name == "":
            raise RuntimeError(f"Feature {self} attribute '{attr_name}' is empty.\n")

        return f"{self.seqid}\t{self.start}\t{self.end}\t{ft_name}\t{self.strand}"

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
        return [
            feature for feature in self.features if feature.featuretype == feature_type
        ]

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
        help="Path to the input GFF file containing genomic features",
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
        dest="p_s_buffer_features",
        type=int,
        default=0,
        required=False,
        help="The number of base pairs to add to the start and end of each feature",
    )
    parser.add_argument(
        "--adjacent_buffer_length",
        dest="p_s_adjacent_buffer_length",
        type=int,
        default=0,
        required=False,
        help="The number of base pairs separating adjacent features between features transcribed by RNA polymerase II and III",
    )
    parser.add_argument(
        "--sheet_name",
        dest="p_i_sheet_name",
        type=str,
        required=True,
        help="The name of the sheet in xlsx_in containing feature information",
    )
    parser.add_argument(
        "--nucleosome_data",
        dest="p_i_nucleosome_data",
        type=str,
        required=False,
        help="The path to the nucleosome data file",
    )
    parser.add_argument(
        "--drop_unmatched",
        dest="p_s_drop_unmatched",
        action="store_true",
        default=False,
        help="Drop features not matched in the nucleosome data from the output",
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
    strand_map = {1: "+", -1: "-", "positive": "+", "negative": "-"}
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

        if feature.strand in strand_map:
            feature.strand = strand_map.get(feature.strand)
        try:
            target_features[row["RNA.Polymerase"]].append(feature)
        except KeyError:
            target_features[row["RNA.Polymerase"]] = [feature]
    return target_features


def build_chrom_id_index(gff_db: FeatureDB) -> Dict[str, str]:
    chrom_id_map = {}
    chrom_prefix = "chr"
    for ft_inst in gff_db.all_features(featuretype="region"):
        try:
            chrom = ft_inst.attributes.get("Name")[0]
            chrom_id_map[ft_inst.chrom] = chrom_prefix + chrom
        except KeyError:
            raise ValueError(
                f"Region feature must have a 'Name' attribute.\n{print(ft_inst)}"
            )
    return chrom_id_map


def load_nucleosome_data(xlsx_file: str, **kwargs) -> Dict[str, Feature]:

    df = pd.read_excel(xlsx_file, sheet_name=kwargs.get("sheet_name"))

    if isinstance(df, dict):
        df = pd.concat(df.values(), ignore_index=True)

    feature_dict = {}
    for index, row in df.iterrows():
        feature = Feature(
            featuretype="ORF",
            seqid=row["Chr"],
            strand=row["Strand"],
            start=int(row["TSS"]),
            end=int(row["TTS"]),
        )
        feature_dict[row["ORF"]] = feature
    return feature_dict


def replace_feature_coordinates(
    nucleosome_data: Dict[str, Feature], features_in: List[Feature]
) -> List[Feature]:
    features_out = []
    for feature in features_in:
        if not feature.id:
            continue

        try:
            orf = re.sub(pattern=r"^(gene|rna)-", repl="", string=feature.id)
            nucleosome_feature = nucleosome_data[orf]
            feature.start = nucleosome_feature.start
            feature.end = nucleosome_feature.end
            if feature.strand == "-":
                feature.start, feature.end = (
                    nucleosome_feature.end,
                    nucleosome_feature.start,
                )
            features_out.append(feature)
        except KeyError:
            # if acc > 0:
            #     raise RuntimeWarning(
            #         f"Feature {feature.id} not found in nucleosome data.\n{print(feature)}"
            #     )
            pass

    print(
        f"Replaced {len(features_out)}/{len(features_in)} feature coordinates with nucleosome data"
    )
    return features_out


def match_features_in_db(
    feature_db: FeatureDB | gff_FeatureDB, features_in: List[Feature]
) -> List[Feature]:
    features_out = []
    chrom_id_map = build_chrom_id_index(feature_db)
    for feature in features_in:
        if feature.has_coordinates():
            features_out.append(feature)
        else:
            for db_feat in feature_db.all_features(featuretype=feature.featuretype):
                db_feat.seqid = chrom_id_map.get(db_feat.seqid)
                if len(feature.attributes):
                    for attr_name, attr_value in feature.attributes.items():
                        if db_feat.attributes.get(attr_name)[0] == attr_value:
                            features_out.append(db_feat)
                else:
                    features_out.append(db_feat)
    return features_out


def load_buffer_specs(
    p_i_xlsx: str, p_i_sheet_name: str, **kwargs
) -> Dict[str, BUFFER_SPEC]:
    buffer_specs = {}
    df = pd.read_excel(p_i_xlsx, sheet_name=p_i_sheet_name)
    for _index, row in df.iterrows():
        if pd.isna(row["Buffer.length"]):
            bl = kwargs.get("default_buffer_length_bp", 0)
        else:
            bl = row["Buffer.length"]

        buffer_specs[row["Feature.name"]] = BUFFER_SPEC(
            feature_name=row["Feature.name"], prefix="", sep=" ", length=bl
        )
    return buffer_specs


def buffer_features(
    features_in: List[Feature], buffer_specs: Dict[str, BUFFER_SPEC], **kwargs
) -> List[Feature]:

    buffer_prefix = kwargs.get("buffer_prefix", ["upstream", "downstream"])
    buffer_sep = kwargs.get("buffer_sep", " ")
    buffer_feature_attr = kwargs.get("buffer_feature_attr", "featuretype")

    features_out = []
    for feature in features_in:
        buffer_spec = buffer_specs.get(feature.featuretype)
        if buffer_spec is None:
            features_out.append(feature)
            continue

        try:
            attr = feature.__getattribute__(buffer_feature_attr)
        except AttributeError as e:
            raise ValueError(
                f"Feature {feature} does not have the attribute {buffer_feature_attr}\n{e}"
            )
        us_feat = Feature(
            id=f"{buffer_prefix[0]}{buffer_sep}{attr}",
            featuretype=attr,
            seqid=feature.seqid,
            start=max(feature.start - buffer_spec.length, 1),
            end=feature.start,
            strand=feature.strand,
        )
        ds_feat = Feature(
            id=f"{buffer_prefix[1]}{buffer_sep}{attr}",
            featuretype=attr,
            seqid=feature.seqid,
            start=feature.end,
            end=feature.end + buffer_spec.length,
            strand=feature.strand,
        )
        features_out.extend([us_feat, feature, ds_feat])
    return features_out


def write_features_to_bed_file(
    p_o_bed: str, features_in: List[Feature], **kwargs
) -> None:
    chrom_prefix = kwargs.get("chrom_prefix", "chr")
    chrom_skip = kwargs.get("chrom_skip", ["MT"])
    chrom_skip.extend([chrom_prefix + chrom for chrom in chrom_skip])
    attr_name = kwargs.get("attr_name", "featuretype")
    with open(p_o_bed, "w") as f:
        for feature in features_in:
            if feature.seqid in chrom_skip:
                continue

            try:
                f.write(f"{feature.to_bed(attr_name)}\n")
            except AttributeError:
                f.write(f"{format_feature_as_bed(feature)}\n")
    return


def overlap_features(feature_a: Feature, feature_b: Feature, **kwargs) -> bool:
    buffer_length = kwargs.get("buffer_length", 0)
    if feature_a.seqid != feature_b.seqid:
        return False
    if (
        feature_a.start - buffer_length > feature_b.end
        or feature_b.start - buffer_length > feature_a.end
    ):
        return False
    return True


def convert_chromosome_id_to_integer(seqid: str) -> int:
    roman_to_int = {"I": 1, "V": 5, "X": 10, "L": 50, "C": 100, "D": 500, "M": 1000}
    roman_id = seqid.replace("chr", "")
    integer_id = 0
    prev_value = 0

    for char in reversed(roman_id):
        current_value = roman_to_int.get(char, 0)
        if current_value == 0:
            if char.isdigit():
                try:
                    return int(roman_id)
                except ValueError:
                    raise ValueError(f"Invalid chromosome ID: {seqid}")
            else:
                raise ValueError(f"Invalid chromosome ID: {seqid}")

        if current_value < prev_value:
            integer_id -= current_value
        else:
            integer_id += current_value
        prev_value = current_value

    return integer_id


def filter_adjacent_features(
    poison_features: List[Feature], query_features: List[Feature], **kwargs
) -> List[Feature]:
    buffer_length = kwargs.get("buffer_length", 0)
    skip_chrom = kwargs.get("skip_chrom", ["MT"])

    skip_chrom.extend([f"chr{chrom}" for chrom in skip_chrom])
    poison_features = [feat for feat in poison_features if feat.seqid not in skip_chrom]
    query_features = [feat for feat in query_features if feat.seqid not in skip_chrom]
    if len(poison_features) == 0 or len(query_features) == 0:
        raise RuntimeError(
            f"Poison features or query features are empty.\n"
            f"Poison features: {len(poison_features)}\n"
            f"Query features: {len(query_features)}"
        )

    # Sort both lists by seqid and start position
    poison_features = sorted(poison_features, key=lambda x: (x.seqid, x.start))
    query_features = sorted(query_features, key=lambda x: (x.seqid, x.start))

    # Initialize indices for both lists
    i = 0
    j = 0

    # Iterate through both lists until one is exhausted
    while i < len(query_features) and j < len(poison_features):
        feat_q = query_features[i]
        feat_p = poison_features[j]
        seqid_q = convert_chromosome_id_to_integer(feat_q.seqid)
        seqid_p = convert_chromosome_id_to_integer(feat_p.seqid)
        if seqid_q < seqid_p:
            i += 1
            continue
        if seqid_q > seqid_p:
            j += 1
            continue

        # If the chromosome IDs are equal, compare the start positions
        if seqid_q == 2 or seqid_p == 2:
            pass
        if feat_p.start - buffer_length > feat_q.end:
            i += 1
        elif feat_q.start - buffer_length > feat_p.end:
            j += 1
        elif overlap_features(feat_q, feat_p, buffer_length=buffer_length):
            query_features.pop(i)
        else:
            i += 1
            j += 1

    return query_features


def main(args=None):
    opts = get_options(args)
    if opts.p_i_gff:
        feature_db = create_gff_db(opts.p_i_gff)
    elif opts.p_i_bed:
        feature_db = load_features_from_bed_file(opts.p_i_bed)
    else:
        raise ValueError("No input file provided")

    target_features = load_target_feature_info(opts.p_i_xlsx, opts.p_i_sheet_name)
    buffer_specs = load_buffer_specs(
        opts.p_i_xlsx, opts.p_i_sheet_name, default_buffer_length_bp=1000
    )

    prev_features = []
    for rna_pol, settings in RNA_POL_SETTINGS.items():
        feature_set = target_features.get(rna_pol)
        if feature_set is None:
            raise RuntimeError(
                f"RNA polymerase {rna_pol} not found in the input file.\n"
                f"Please check the input file and try again."
            )

        output_bed_path = get_bed_file_name(opts.p_o_bed_prefix, rna_pol)
        feature_set = match_features_in_db(feature_db, feature_set)

        if rna_pol == 3:
            prev_features.extend(feature_set)
            feature_set = buffer_features(feature_set, buffer_specs)
            # Fix feature IDs for RNA Pol III features specifically
            for feature in feature_set:
                if feature.id == "":
                    feature.id = feature.featuretype
                elif re.match(pattern=r"^(gene|rna)-", string=feature.id):
                    feature.id = feature.featuretype
                else:
                    pass
        else:
            if opts.p_i_nucleosome_data:
                nucleosome_data = load_nucleosome_data(opts.p_i_nucleosome_data)
                nucl_features = replace_feature_coordinates(
                    nucleosome_data, feature_set
                )
                if opts.p_s_drop_unmatched:
                    feature_set = nucl_features

            feature_set = buffer_features(
                feature_set,
                buffer_specs,
                buffer_feature_attr="id",
                buffer_prefix=["", ""],
                buffer_sep="",
            )
            pre_count = len(feature_set)
            if pre_count == 0:
                raise RuntimeError(f"No features found for RNA Pol {rna_pol}.")

            feature_set = filter_adjacent_features(
                prev_features,
                feature_set,
                buffer_length=opts.p_s_adjacent_buffer_length,
            )
            post_count = len(feature_set)
            print(f"Filtered {pre_count - post_count} adjacent features")

        write_features_to_bed_file(
            p_o_bed=output_bed_path,
            features_in=feature_set,
            attr_name=settings["feature_attribute"],
        )
    return


if __name__ == "__main__":
    main()
