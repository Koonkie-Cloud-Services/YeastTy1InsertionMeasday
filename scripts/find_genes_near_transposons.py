#!/usr/bin/env python

"""Identify features within a specified distance of Ty1 elements in a BED file."""

import argparse
import sys
from typing import List, Tuple


def get_options(args=sys.argv[1:]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find genes (non-RNA) within a specified distance of Ty1 elements"
    )
    parser.add_argument(
        "--bed",
        dest="bed_file",
        type=str,
        required=True,
        help="Input BED file with feature locations"
    )
    parser.add_argument(
        "--features",
        dest="features",
        type=str,
        default="gene",
        help="Comma-separated list of feature name patterns to scan (default: 'gene')"
    )
    parser.add_argument(
        "--distance",
        dest="distance",
        type=int,
        default=1000,
        help="Maximum distance in bp to consider as 'near' a Ty1 element (default: 1000)"
    )
    parser.add_argument(
        "--output",
        dest="output_file",
        type=str,
        required=True,
        help="Output file path for results"
    )
    return parser.parse_args(args=args)


def parse_bed_file(bed_file: str) -> List[Tuple[str, int, int, str, str]]:
    """Parse a BED file into a list of (chrom, start, end, name, strand) tuples."""
    records = []
    with open(bed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3]
            strand = fields[4] if len(fields) > 4 else "."
            records.append((chrom, start, end, name, strand))
    return records


def is_ty1_element(name: str) -> bool:
    """Check if a feature name indicates a Ty1 element (retrotransposon or LTR)."""
    return name.startswith("Ty1(")


def matches_feature_patterns(name: str, patterns: List[str]) -> bool:
    """Check if a feature name matches any of the specified patterns.

    A record matches if its name contains any pattern as a substring.
    RNA features (tRNA, snoRNA, snRNA, rRNA, other_RNA) are excluded.
    """
    rna_types = {"tRNA", "snoRNA", "snRNA", "rRNA", "other_RNA"}
    if name in rna_types:
        return False
    # Also exclude names ending with _upstream or _downstream
    if name.endswith("_upstream") or name.endswith("_downstream"):
        return False
    return any(pattern in name for pattern in patterns)


def distance_between(start1: int, end1: int, start2: int, end2: int) -> int:
    """Calculate the distance between two intervals. Returns 0 if they overlap."""
    if start1 >= end2:
        return start1 - end2
    elif start2 >= end1:
        return start2 - end1
    else:
        return 0


def find_features_near_ty1(records: List[Tuple[str, int, int, str, str]],
                           feature_patterns: List[str],
                           max_distance: int) -> List[Tuple[str, int, int, str, str, str, int]]:
    """Find features matching the patterns that are within max_distance of a Ty1 element.

    Returns list of (chrom, start, end, name, strand, nearest_ty1, distance) tuples.
    """
    # Index Ty1 elements by chromosome
    ty1_by_chrom = {}
    for chrom, start, end, name, strand in records:
        if is_ty1_element(name):
            ty1_by_chrom.setdefault(chrom, []).append((start, end, name))

    # Find matching features near Ty1 elements
    results = []
    seen = set()
    for chrom, start, end, name, strand in records:
        if not matches_feature_patterns(name, feature_patterns):
            continue
        if chrom not in ty1_by_chrom:
            continue

        for ty1_start, ty1_end, ty1_name in ty1_by_chrom[chrom]:
            dist = distance_between(start, end, ty1_start, ty1_end)
            if dist <= max_distance:
                key = (chrom, start, end, name, strand)
                if key not in seen:
                    seen.add(key)
                    results.append((chrom, start, end, name, strand, ty1_name, dist))
                break

    return results


def write_results(results: List[Tuple[str, int, int, str, str, str, int]],
                  output_file: str) -> None:
    """Write results to a tab-separated file."""
    with open(output_file, "w") as fh:
        fh.write("#chrom\tstart\tend\tfeature\tstrand\tnearest_ty1\tdistance\n")
        for chrom, start, end, name, strand, ty1_name, dist in results:
            fh.write(f"{chrom}\t{start}\t{end}\t{name}\t{strand}\t{ty1_name}\t{dist}\n")


def main(args=None):
    opts = get_options(args)
    feature_patterns = [p.strip() for p in opts.features.split(",")]

    records = parse_bed_file(opts.bed_file)
    results = find_features_near_ty1(records, feature_patterns, opts.distance)

    write_results(results, opts.output_file)
    print(f"Found {len(results)} features within {opts.distance} bp of Ty1 elements")
    print(f"Results written to {opts.output_file}")


if __name__ == "__main__":
    main()
