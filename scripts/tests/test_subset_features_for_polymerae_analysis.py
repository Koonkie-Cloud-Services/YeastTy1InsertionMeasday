import os
import unittest

class Tester(unittest.TestCase):
    def setUp(self):
        self.test_bed_out_prefix = "tests/test_out"
        return
    
    def tearDown(self):
        from subset_features_for_polymerase_analysis import get_bed_file_name
        for i in range(1, 4):
            test_bed = get_bed_file_name(self.test_bed_out_prefix, i)
            if os.path.exists(test_bed):
                os.remove(test_bed)
        return super().tearDown()
    
    def test_main(self):
        from subset_features_for_polymerase_analysis import main
        main([
        "--bed_in",
        "../data/output/feature_locations/Saccharomyces_cerevisiae_S288C.bed",
        "--xlsx_in",
        "../data/input/custom_feature_annotations.xlsx",
        "--sheet_name",
        "polymerase_analysis",
        "--bed_out_prefix",
        self.test_bed_out_prefix])
        return
    