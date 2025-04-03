import os
import unittest

class Tester(unittest.TestCase):
    def setUp(self):
        self.test_bed_out_prefix = "tests/results/test_out"
        return
    
    # def test_subset_positions(self):
    #     from intersect_coverage_and_feature_beds import subset_positions
    #     import pandas as pd

    #     with self.subTest("Subset positions from BED files"):
    #         cov_df = pd.DataFrame({
    #             'seqid': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
    #             'strand': ['+', '+', '+', '+', '+'],
    #             'position': [1, 2, 3, 4, 5],
    #             'coverage': [0, 0, 1, 4, 2]
    #         })
    #         feat_df = pd.DataFrame({
    #             'seqid': ['chr1', 'chr1', 'chr1'],
    #             'start': [2, 3, 5],
    #             'end': [2, 4, 5],
    #             'featuretype': ['upstream gene', 'gene', 'downstream gene'],
    #             'strand': ['+', '+', '+']
    #         })
    #         res_df = subset_positions(cov_df, feat_df)
    #         self.assertEqual(res_df.shape[0], 4)
    #         self.assertEqual(res_df.shape[1], 5)
    #         self.assertEqual(res_df['position'].tolist(), [2, 3, 4, 5])
    #         self.assertEqual(res_df['featuretype'].tolist(), ['upstream gene', 'gene', 'gene', 'downstream gene'])
    #     return
    
    def test_intersect_bed_files(self):
        from intersect_coverage_and_feature_beds import intersect_bed_files

        with self.subTest("Subsetting features from GFF file"):
            intersect_bed_files([
            "--feature_bed",
            "tests/data/genes.bed",
            "--coverage_bed",
            "tests/data/coverage.bed",
            "--output_bed",
            f"{self.test_bed_out_prefix}.bed"])
        return