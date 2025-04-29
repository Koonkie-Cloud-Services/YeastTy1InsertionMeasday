import os
import unittest

from subset_features_for_polymerase_analysis import Feature

class TestMain(unittest.TestCase):
    def setUp(self):
        self.test_bed_out_prefix = "tests/results/TestMain"
        self.nucelosome_xlsx = "../data/input/13059_2018_1398_MOESM2_ESM.xlsx"
        if not os.path.exists(os.path.dirname(self.test_bed_out_prefix)):
            os.makedirs(os.path.dirname(self.test_bed_out_prefix))
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

        with self.subTest("Subsetting features from GFF file"):
            main([
            "--gff_in",
            "../data/output/feature_locations/Saccharomyces_cerevisiae_S288C.gff",
            "--xlsx_in",
            "../data/input/custom_feature_annotations.xlsx",
            "--sheet_name",
            "polymerase_analysis",
            "--bed_out_prefix",
            self.test_bed_out_prefix,
            "--nucleosome_data",
            self.nucelosome_xlsx])

        # with self.subTest("Subsetting features from BED file"):
        #     main([
        #     "--bed_in",
        #     "../data/output/feature_locations/Saccharomyces_cerevisiae_S288C.bed",
        #     "--xlsx_in",
        #     "../data/input/custom_feature_annotations.xlsx",
        #     "--sheet_name",
        #     "polymerase_analysis",
        #     "--bed_out_prefix",
        #     self.test_bed_out_prefix])
        return
    
class TestNucleosomeIntegration(unittest.TestCase):
    def setUp(self):
        self.test_dir = "tests/results/TestNucleosomeIntegration"
        self.nucelosome_xlsx = "../data/input/13059_2018_1398_MOESM2_ESM.xlsx"
    
    def test_load_nucleosome_data(self):
        from subset_features_for_polymerase_analysis import load_nucleosome_data

        result = load_nucleosome_data(self.nucelosome_xlsx)

        self.assertIsInstance(result, dict)
        self.assertEqual(len(result), 5542)
        self.assertEqual(result["YAL068C"].start, 2169)
        return

    def test_replace_nucleosome_data(self):
        from subset_features_for_polymerase_analysis import replace_feature_coordinates

        nucleosome_data = {
            "YAL068C": Feature(featuretype="nucleosome", seqid="chrI", start=2169, end=2269, strand="+"),
            "YAL069W": Feature(featuretype="nucleosome", seqid="chrI", start=2269, end=2369, strand="+"),
        }
        feature = Feature(featuretype="gene", seqid="chrI", start=2000, end=3000, strand="+", id="gene-YAL068C")
        replace_feature_coordinates(nucleosome_data, [feature])

        self.assertEqual(feature.start, 2169)
        self.assertEqual(feature.end, 2269)
        return

class TestFilterAdjacentFeatures(unittest.TestCase):
    def setUp(self):
        self.poison_features = [
            Feature(featuretype="gene", seqid="chrI", start=100, end=200, strand="+"),
            Feature(featuretype="gene", seqid="chrI", start=300, end=400, strand="+"),
            Feature(featuretype="gene", seqid="chrII", start=100, end=200, strand="+"),
        ]
        self.query_features = [
            Feature(featuretype="gene", seqid="chrI", start=150, end=250, strand="+"),
            Feature(featuretype="gene", seqid="chrI", start=500, end=600, strand="+"),
            Feature(featuretype="gene", seqid="chrII", start=50, end=150, strand="+"),
        ]

    def test_filter_adjacent_features_iteration(self):
        from subset_features_for_polymerase_analysis import filter_adjacent_features
        buffer_length = 50
        filtered_features = filter_adjacent_features(
            self.poison_features, self.query_features, buffer_length=buffer_length
        )

        # Ensure the iteration correctly filters overlapping features
        self.assertEqual(len(filtered_features), 1)
        self.assertEqual(filtered_features[0].start, 500)
        self.assertEqual(filtered_features[0].end, 600)

    def test_convert_chromosome_id_to_integer(self):
        from subset_features_for_polymerase_analysis import convert_chromosome_id_to_integer

        with self.subTest("Valid chromosome ID"):
            chromosome_id = "chrI"
            expected_result = 1
            result = convert_chromosome_id_to_integer(chromosome_id)
            self.assertEqual(result, expected_result)

        with self.subTest("Invalid chromosome ID"):
            chromosome_id = "chrY"
            self.assertRaises(ValueError,
                              lambda: convert_chromosome_id_to_integer(chromosome_id))

        with self.subTest("Chromosome ID with arabic number"):
            chromosome_id = "chr3"
            expected_result = 3
            result = convert_chromosome_id_to_integer(chromosome_id)
            self.assertEqual(result, expected_result)
        
        return

if __name__ == "__main__":
    unittest.main()