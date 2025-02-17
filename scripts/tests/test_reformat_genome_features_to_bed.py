import unittest

import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reformat_genome_features_to_bed import parse_gff_to_bed, create_gff_db


class TestParseGFFToBED(unittest.TestCase):

    def setUp(self):
        self.gff_file = "test_input.gff"
        self.bed_file = "test_output.bed"
        self.gff_content = """##gff-version 3
NC_001133.9\tRefSeq\tregion\t1\t230218\t.\t+\t.\tID=NC_001133.9:1..230218;Dbxref=taxon:559292;Name=I;chromosome=I;gbkey=Src;genome=chromosome;mol_type=genomic DNA;strain=S288C
NC_001133.9\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=Gene1
NC_001133.9\t.\ttRNA\t3000\t4000\t.\t-\t.\tID=tRNA1;Name=tRNA1
NC_001133.9\t.\trepeat_region\t5000\t6000\t.\t+\t.\tID=repeat1;Name=Repeat1
NC_001133.9\tRefSeq\tcentromere\t151465\t151582\t.\t+\t.\tID=id-NC_001133.9:151465..151582;Dbxref=SGD:S000028887;Name=CEN1;chromosome=I
NC_001133.9\tRefSeq\tcentromere\t151465\t151474\t.\t+\t.\tID=id-NC_001133.9:151465..151474;Note=CEN1
NC_001133.9\tRefSeq\tcentromere\t151475\t151557\t.\t+\t.\tID=id-NC_001133.9:151475..151557;Note=CEN1
NC_001133.9\tRefSeq\tcentromere\t151558\t151582\t.\t+\t.\tID=id-NC_001133.9:151558..151582;Note=CEN1
NC_001224.1\tRefSeq\tregion\t1\t85779\t.\t+\t.\tID=NC_001224.1:1..85779;Dbxref=taxon:559292;Name=MT;gbkey=Src;genome=mitochondrion;mol_type=genomic DNA;strain=S288C;Is_circular=true;substrain=FY1679
NC_001224.1\tRefSeq\tgene\t731\t802\t.\t+\t.\tID=gene-tP(UGG)Q;Dbxref=GeneID:854578,SGD:S000007307;Name=tP(UGG)Q;gbkey=Gene;gene_biotype=tRNA
        """
        self.expected_bed_content = """chrI\t1000\t2000\tgene1(gene)\t+
chrI\t3000\t4000\ttRNA\t-
chrI\t5000\t6000\tLTR\t+
"""
        with open(self.gff_file, 'w') as fh:
            fh.write(self.gff_content.strip())
        self.gff_db = create_gff_db(gff_file = self.gff_file)
        assert self.gff_db.count_features_of_type() == 10
        return

    def tearDown(self):
        for test_file in [self.gff_file, self.bed_file]:
            if os.path.exists(test_file):
                os.remove(test_file)
        return super().tearDown()

    def test_parse_gff_to_bed(self):
        parse_gff_to_bed(self.gff_db,
                         self.bed_file,
                         feature_types_to_skip=["centromere", "region"])

        with open(self.bed_file, "r") as f:
            result = f.read()

        self.assertEqual(result, self.expected_bed_content)
        return

    def test_remove_redundant_features(self):
        from reformat_genome_features_to_bed import remove_redundant_features
        remove_redundant_features(self.gff_db)
        self.assertEqual(7, self.gff_db.count_features_of_type())
        self.assertEqual(1, self.gff_db.count_features_of_type("centromere"))
        return

if __name__ == "__main__":
    unittest.main()
