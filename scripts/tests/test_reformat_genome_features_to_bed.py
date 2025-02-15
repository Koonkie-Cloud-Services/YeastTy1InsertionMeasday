import unittest

import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reformat_genome_features_to_bed import parse_gff_to_bed, create_gff_db
import gffutils


class TestParseGFFToBED(unittest.TestCase):

    def setUp(self):
        self.gff_file = "test_input.gff"
        self.bed_file = "test_output.bed"
        self.gff_content = """##gff-version 3
NC_001133.9\tRefSeq\tregion\t1\t230218\t.\t+\t.\tID=NC_001133.9:1..230218;Dbxref=taxon:559292;Name=I;chromosome=I;gbkey=Src;genome=chromosome;mol_type=genomic DNA;strain=S288C
chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=Gene1
chr1\t.\ttRNA\t3000\t4000\t.\t-\t.\tID=tRNA1;Name=tRNA1
chr1\t.\trepeat_region\t5000\t6000\t.\t+\t.\tID=repeat1;Name=Repeat1
    """
        with open(self.gff_file, 'w') as fh:
            fh.write(self.gff_content)
        db_fn = create_gff_db(gff_file = self.gff_file)
        self.gff_db = gffutils.FeatureDB(db_fn)
        self.expected_bed_content = """1000\t2000\tgene1(gene)\t+
3000\t4000\ttRNA\t-
5000\t6000\tLTR\t+
"""

    def tearDown(self):
        for test_file in [self.gff_file, self.bed_file]:
            os.remove(test_file)
        return super().tearDown()

    def test_parse_gff_to_bed(self):
        parse_gff_to_bed(self.gff_db, self.bed_file)

        with open(self.bed_file, "r") as f:
            result = f.read()

        self.assertEqual(result, self.expected_bed_content)


if __name__ == "__main__":
    unittest.main()
