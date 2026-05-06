#!/usr/bin/env python3
"""
Unit tests for calculate_fractional_coverage.py

Tests fractional read assignment logic with various scenarios:
- Equal split across multiple features
- Strand-specific matching
- Single overlaps
- Edge cases
"""

import unittest
import sys
import io
from io import StringIO


class TestFractionalCoverage(unittest.TestCase):
    """Test cases for fractional coverage calculation."""
    
    def setUp(self):
        """Import the module to test."""
        # Add parent directory to path to import the script
        import os
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        if parent_dir not in sys.path:
            sys.path.insert(0, parent_dir)
        import calculate_fractional_coverage
        self.module = calculate_fractional_coverage
    
    def test_equal_split_three_features(self):
        """Test that one read overlapping 3 features gives 0.333... to each."""
        # Input: one read (read1) overlapping 3 features on same strand
        input_data = """chrI\t100\t150\tread1\t255\t+\tchrI\t90\t160\tgene1\t+
chrI\t100\t150\tread1\t255\t+\tchrI\t120\t180\tgene2\t+
chrI\t100\t150\tread1\t255\t+\tchrI\t140\t200\tgene3\t+"""
        
        input_stream = StringIO(input_data)
        feature_counts = self.module.calculate_fractional_coverage(input_stream)
        
        # Each feature should get 1/3 of the read
        self.assertEqual(len(feature_counts), 3)
        self.assertAlmostEqual(feature_counts[('chrI', '90', '160', 'gene1', '+')], 0.3333, places=4)
        self.assertAlmostEqual(feature_counts[('chrI', '120', '180', 'gene2', '+')], 0.3333, places=4)
        self.assertAlmostEqual(feature_counts[('chrI', '140', '200', 'gene3', '+')], 0.3333, places=4)
        
        # Total should sum to 1
        total = sum(feature_counts.values())
        self.assertAlmostEqual(total, 1.0, places=6)
    
    def test_single_overlap(self):
        """Test that one read overlapping one feature gives count of 1.0."""
        input_data = """chrI\t100\t150\tread1\t255\t+\tchrI\t90\t160\tgene1\t+"""
        
        input_stream = StringIO(input_data)
        feature_counts = self.module.calculate_fractional_coverage(input_stream)
        
        self.assertEqual(len(feature_counts), 1)
        self.assertEqual(feature_counts[('chrI', '90', '160', 'gene1', '+')], 1.0)
    
    def test_multiple_reads_same_feature(self):
        """Test that multiple reads to the same feature sum correctly."""
        # Two reads, each uniquely overlapping gene1
        input_data = """chrI\t100\t150\tread1\t255\t+\tchrI\t90\t160\tgene1\t+
chrI\t200\t250\tread2\t255\t+\tchrI\t90\t160\tgene1\t+"""
        
        input_stream = StringIO(input_data)
        feature_counts = self.module.calculate_fractional_coverage(input_stream)
        
        self.assertEqual(len(feature_counts), 1)
        self.assertEqual(feature_counts[('chrI', '90', '160', 'gene1', '+')], 2.0)
    
    def test_mixed_overlap_patterns(self):
        """Test complex scenario with different overlap patterns."""
        # read1: overlaps gene1 only (1.0)
        # read2: overlaps gene2 and gene3 (0.5 each)
        # read3: overlaps gene3 only (1.0)
        # Expected: gene1=1.0, gene2=0.5, gene3=1.5
        input_data = """chrI\t100\t150\tread1\t255\t+\tchrI\t90\t160\tgene1\t+
chrI\t200\t250\tread2\t255\t+\tchrI\t190\t260\tgene2\t+
chrI\t200\t250\tread2\t255\t+\tchrI\t240\t300\tgene3\t+
chrI\t300\t350\tread3\t255\t+\tchrI\t240\t300\tgene3\t+"""
        
        input_stream = StringIO(input_data)
        feature_counts = self.module.calculate_fractional_coverage(input_stream)
        
        self.assertEqual(len(feature_counts), 3)
        self.assertEqual(feature_counts[('chrI', '90', '160', 'gene1', '+')], 1.0)
        self.assertEqual(feature_counts[('chrI', '190', '260', 'gene2', '+')], 0.5)
        self.assertEqual(feature_counts[('chrI', '240', '300', 'gene3', '+')], 1.5)
        
        # Total reads should sum to 3
        total = sum(feature_counts.values())
        self.assertAlmostEqual(total, 3.0, places=6)
    
    def test_strand_separation(self):
        """Test that reads only count toward features on matching strand."""
        # read1 (+) overlaps gene1 (+) -> counts (1.0)
        # read2 (-) overlaps gene2 (-) -> counts (1.0)
        # read3 (+) overlaps gene3 (-) -> filtered out (mismatched strand)
        input_data = """chrI\t100\t150\tread1\t255\t+\tchrI\t90\t160\tgene1\t+
chrI\t200\t250\tread2\t255\t-\tchrI\t190\t260\tgene2\t-
chrI\t300\t350\tread3\t255\t+\tchrI\t290\t360\tgene3\t-"""
        
        input_stream = StringIO(input_data)
        feature_counts = self.module.calculate_fractional_coverage(input_stream)
        
        # Only 2 features should have counts (strand-matched ones)
        self.assertEqual(len(feature_counts), 2)
        self.assertEqual(feature_counts[('chrI', '90', '160', 'gene1', '+')], 1.0)
        self.assertEqual(feature_counts[('chrI', '190', '260', 'gene2', '-')], 1.0)
        # gene3 should not be in the results (read3 had wrong strand)
        self.assertNotIn(('chrI', '290', '360', 'gene3', '-'), feature_counts)
    
    def test_strand_mismatch_filtering(self):
        """Test that strand mismatches are properly filtered."""
        # read1 overlaps 2 features: gene1 on + (matches), gene2 on - (doesn't match)
        # read1 should only get counted toward gene1
        input_data = """chrI\t100\t150\tread1\t255\t+\tchrI\t90\t160\tgene1\t+
chrI\t100\t150\tread1\t255\t+\tchrI\t120\t180\tgene2\t-"""
        
        input_stream = StringIO(input_data)
        feature_counts = self.module.calculate_fractional_coverage(input_stream)
        
        # Only gene1 should receive the full count (1.0)
        self.assertEqual(len(feature_counts), 1)
        self.assertEqual(feature_counts[('chrI', '90', '160', 'gene1', '+')], 1.0)
        # gene2 should not be in results
        self.assertNotIn(('chrI', '120', '180', 'gene2', '-'), feature_counts)
    
    def test_empty_input(self):
        """Test that empty input produces empty output."""
        input_data = ""
        
        input_stream = StringIO(input_data)
        feature_counts = self.module.calculate_fractional_coverage(input_stream)
        
        self.assertEqual(len(feature_counts), 0)
    
    def test_parse_intersect_line(self):
        """Test parsing of bedtools intersect -wa -wb output."""
        line = "chrI\t100\t150\tread1\t255\t+\tchrI\t90\t160\tgene1\t+"
        
        read_id, read_strand, f_chrom, f_start, f_end, f_id, f_strand = \
            self.module.parse_intersect_line(line)
        
        self.assertEqual(read_id, 'read1')
        self.assertEqual(read_strand, '+')
        self.assertEqual(f_chrom, 'chrI')
        self.assertEqual(f_start, '90')
        self.assertEqual(f_end, '160')
        self.assertEqual(f_id, 'gene1')
        self.assertEqual(f_strand, '+')
    
    def test_malformed_line_handling(self):
        """Test that malformed lines are handled gracefully."""
        input_data = """chrI\t100\t150\tread1\t255\t+\tchrI\t90\t160\tgene1\t+
chrI\t100\t150\tmalformed_line
chrI\t200\t250\tread2\t255\t+\tchrI\t190\t260\tgene2\t+"""
        
        input_stream = StringIO(input_data)
        # Should process valid lines and skip malformed ones
        feature_counts = self.module.calculate_fractional_coverage(input_stream)
        
        self.assertEqual(len(feature_counts), 2)
        self.assertEqual(feature_counts[('chrI', '90', '160', 'gene1', '+')], 1.0)
        self.assertEqual(feature_counts[('chrI', '190', '260', 'gene2', '+')], 1.0)
    
    def test_four_way_overlap(self):
        """Test read overlapping 4 features (each gets 0.25)."""
        input_data = """chrI\t100\t150\tread1\t255\t+\tchrI\t90\t160\tgene1\t+
chrI\t100\t150\tread1\t255\t+\tchrI\t120\t180\tgene2\t+
chrI\t100\t150\tread1\t255\t+\tchrI\t140\t200\tgene3\t+
chrI\t100\t150\tread1\t255\t+\tchrI\t145\t210\tgene4\t+"""
        
        input_stream = StringIO(input_data)
        feature_counts = self.module.calculate_fractional_coverage(input_stream)
        
        self.assertEqual(len(feature_counts), 4)
        for key in feature_counts:
            self.assertEqual(feature_counts[key], 0.25)
        
        # Total should sum to 1
        total = sum(feature_counts.values())
        self.assertAlmostEqual(total, 1.0, places=6)


class TestOutputFormatting(unittest.TestCase):
    """Test output formatting functions."""
    
    def setUp(self):
        """Import the module to test."""
        # Add parent directory to path to import the script
        import os
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        if parent_dir not in sys.path:
            sys.path.insert(0, parent_dir)
        import calculate_fractional_coverage
        self.module = calculate_fractional_coverage
    
    def test_write_output(self):
        """Test that output is formatted correctly."""
        feature_counts = {
            ('chrI', '90', '160', 'gene1', '+'): 1.3333,
            ('chrI', '120', '180', 'gene2', '+'): 0.5000,
            ('chrII', '200', '300', 'gene3', '-'): 2.0000,
        }
        
        # Capture stdout
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()
        
        try:
            self.module.write_output(feature_counts)
            output = sys.stdout.getvalue()
        finally:
            sys.stdout = old_stdout
        
        lines = output.strip().split('\n')
        self.assertEqual(len(lines), 3)
        
        # Check sorting (chrI before chrII, then by start position)
        self.assertTrue(lines[0].startswith('chrI\t90\t160\tgene1\t+\t1.3333'))
        self.assertTrue(lines[1].startswith('chrI\t120\t180\tgene2\t+\t0.5000'))
        self.assertTrue(lines[2].startswith('chrII\t200\t300\tgene3\t-\t2.0000'))


def run_tests():
    """Run all tests."""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    suite.addTests(loader.loadTestsFromTestCase(TestFractionalCoverage))
    suite.addTests(loader.loadTestsFromTestCase(TestOutputFormatting))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_tests()
    sys.exit(0 if success else 1)
