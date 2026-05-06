#!/usr/bin/env python3
"""
Calculate fractional read coverage for genomic features.

This script processes read-feature intersections and assigns fractional counts
when reads overlap multiple features. Each read is split equally among all
overlapping features (e.g., 1 read overlapping 3 features = 0.333... to each).

Input: Output from 'bedtools intersect -wa -wb' (TSV format)
  - Columns 1-6: Read info (chrom, start, end, read_id, score, strand)
  - Columns 7-11: Feature info (chrom, start, end, feature_id, strand)

Output: Feature-level coverage (TSV format)
  - Columns: chrom, start, end, feature_id, strand, fractional_count

Usage:
  bedtools intersect -wa -wb -a reads.bed -b features.bed | \
    python calculate_fractional_coverage.py > output.bed
"""

import sys
from collections import defaultdict
from typing import Dict, List, Tuple


def parse_intersect_line(line: str) -> Tuple[str, str, str, str, str, str, str]:
    """
    Parse a line from bedtools intersect -wa -wb output.
    
    Returns:
        (read_id, read_strand, feature_chrom, feature_start, feature_end, feature_id, feature_strand)
    """
    fields = line.strip().split('\t')
    
    if len(fields) < 11:
        raise ValueError(f"Expected at least 11 columns, got {len(fields)}: {line}")
    
    # Read info (columns 0-5)
    read_id = fields[3]
    read_strand = fields[5]
    
    # Feature info (columns 6-10)
    feature_chrom = fields[6]
    feature_start = fields[7]
    feature_end = fields[8]
    feature_id = fields[9]
    feature_strand = fields[10]
    
    return (read_id, read_strand, feature_chrom, feature_start, feature_end, feature_id, feature_strand)


def calculate_fractional_coverage(intersect_stream) -> Dict[Tuple[str, str, str, str, str], float]:
    """
    Calculate fractional coverage from read-feature intersections.
    
    For each read, assign equal fractional counts (1/N) to all N overlapping features.
    
    Args:
        intersect_stream: Input stream of bedtools intersect output
        
    Returns:
        Dictionary mapping (chrom, start, end, feature_id, strand) -> fractional_count
    """
    # First pass: Group intersections by read_id (only if strands match)
    read_to_features: Dict[str, List[Tuple[str, str, str, str, str]]] = defaultdict(list)
    
    for line in intersect_stream:
        if not line.strip():
            continue
            
        try:
            read_id, read_strand, f_chrom, f_start, f_end, f_id, f_strand = parse_intersect_line(line)
            
            # Only count read-feature pairs where strands match
            if read_strand != f_strand:
                continue
            
            feature_key = (f_chrom, f_start, f_end, f_id, f_strand)
            read_to_features[read_id].append(feature_key)
        except (ValueError, IndexError) as e:
            print(f"Warning: Skipping malformed line: {e}", file=sys.stderr)
            continue
    
    # Second pass: Calculate fractional counts
    feature_counts: Dict[Tuple[str, str, str, str, str], float] = defaultdict(float)
    
    for read_id, features in read_to_features.items():
        # Equal split: each feature gets 1/N
        num_features = len(features)
        fractional_count = 1.0 / num_features
        
        for feature_key in features:
            feature_counts[feature_key] += fractional_count
    
    return feature_counts


def write_output(feature_counts: Dict[Tuple[str, str, str, str, str], float]):
    """
    Write fractional coverage results in BED format.
    
    Output format: chrom, start, end, feature_id, strand, fractional_count
    """
    # Sort by chromosome, then start position
    sorted_features = sorted(feature_counts.items(), 
                            key=lambda x: (x[0][0], int(x[0][1])))
    
    for feature_key, count in sorted_features:
        chrom, start, end, feature_id, strand = feature_key
        # Format count with 4 decimal places for readability
        print(f"{chrom}\t{start}\t{end}\t{feature_id}\t{strand}\t{count:.4f}")


def main():
    """Main entry point."""
    if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help']:
        print(__doc__)
        sys.exit(0)
    
    # Read from stdin
    feature_counts = calculate_fractional_coverage(sys.stdin)
    
    # Write to stdout
    write_output(feature_counts)
    
    # Log statistics to stderr
    total_features = len(feature_counts)
    total_reads = sum(feature_counts.values())
    print(f"Processed {total_reads:.4f} reads across {total_features} features", 
          file=sys.stderr)


if __name__ == '__main__':
    main()
