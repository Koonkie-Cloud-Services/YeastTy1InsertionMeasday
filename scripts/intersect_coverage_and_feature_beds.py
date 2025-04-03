import sys
import pandas as pd
import argparse

def parse_arguments(args):
    parser = argparse.ArgumentParser(description='Subset positions from BED files.')
    parser.add_argument('--coverage_bed', dest="p_i_coverage_bed_path", type=str, help='Path to the three-column genome coverage BED file.')
    parser.add_argument('--feature_bed', dest="p_i_feature_bed_path", type=str, help='Path to the five-column BED file with feature positions.')
    parser.add_argument('--output_bed', dest="p_o_bed", type=str, help='Path to the output BED file.')
    return parser.parse_args(args)

def read_tsv(file_path: str, **kwargs) -> pd.DataFrame:
    return pd.read_csv(file_path, sep='\t', **kwargs)

def write_bed_intersect(cov_df: pd.DataFrame, feat_df: pd.DataFrame, bed_path: str) -> None:
    cov_df.columns = ['seqid', 'strand', 'position', 'coverage']
    feat_df.columns = ['seqid', 'start', 'end', 'featuretype', 'strand']

    bed_out = open(bed_path, 'w')

    def load_feature_index(feat_df):
        feature_idx = {}
        # Index the features in feat_df by chromosome, start and end coords
        for _i, row in feat_df.iterrows():
            if row['seqid'] not in feature_idx:
                feature_idx[row['seqid']] = {
                    '+' : {},
                    '-' : {}
                }

            try:
                feature_idx[row['seqid']][row['strand']][row['start']].append((row['featuretype'], row['end']))
            except KeyError:
                feature_idx[row['seqid']][row['strand']][row['start']] = [(row['featuretype'], row['end'])]
        return feature_idx

    feature_idx = load_feature_index(feat_df)

    feature_queue = []
    # drop_indices = []
    prev_seqid = None
    prev_strand = None
    for i, row in cov_df.iterrows():
        if row['seqid'] != prev_seqid or row['strand'] != prev_strand:
            feature_queue = []
            prev_seqid = row['seqid']
            prev_strand = row['strand']

        if row['seqid'] not in feature_idx:
            raise ValueError(f"Chromosome {row['seqid']} not found in feature BED file.")

        try:
            feature_queue.extend(feature_idx[row['seqid']][row['strand']][row['position']])
        except KeyError:
            pass

        j = 0
        while feature_queue and j < len(feature_queue):
            if row['position'] > feature_queue[j][1]:
                feature_queue.pop(j)
                continue
            bed_out.write(row.to_csv(index=False, header=False, lineterminator="\t") + f"{feature_queue[j][0]}\n")
            j += 1

        # if not feature_queue:
        #     drop_indices.append(i)
    
    # cov_df.drop(drop_indices, inplace=True)
    bed_out.close()
    
    return

def intersect_bed_files(args = None):
    if args is None:
        args = sys.argv[1:]
    opts = parse_arguments(args)

    cov_df = read_tsv(opts.p_i_coverage_bed_path, header=None)
    feat_df = read_tsv(opts.p_i_feature_bed_path, header=None)

    write_bed_intersect(cov_df, feat_df, opts.p_o_bed)
    return

if __name__ == "__main__":
    intersect_bed_files()