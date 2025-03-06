import unittest

import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reformat_genome_features_to_bed import create_gff_db


class TestParseGFFToBED(unittest.TestCase):

    def setUp(self):
        self.gff_file = "test_input.gff"
        self.bed_file = "test_output.bed"
        self.gff_content = """##gff-version 3
NC_001133.9\tRefSeq\tregion\t1\t230218\t.\t+\t.\tID=NC_001133.9:1..230218;Dbxref=taxon:559292;Name=I;chromosome=I;gbkey=Src;genome=chromosome;mol_type=genomic DNA;strain=S288C
NC_001133.9\tMock\ttRNA\t640\t700\t.\t+\t.\tID=tRNA1;Name=tRNA1
NC_001133.9\tMock\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=Gene1
NC_001133.9\tMock\ttRNA\t3000\t4000\t.\t-\t.\tID=tRNA2;Name=tRNA2
NC_001133.9\tMock\tlong_terminal_repeat\t5000\t6000\t.\t+\t.\tID=repeat1;Name=Repeat1;Note=Ty1 LTR;gbkey=repeat_region;rpt_type=long_terminal_repeat
NC_001133.9\tRefSeq\tcentromere\t151465\t151582\t.\t+\t.\tID=id-NC_001133.9:151465..151582;Dbxref=SGD:S000028887;Name=CEN1;chromosome=I
NC_001133.9\tRefSeq\tcentromere\t151465\t151474\t.\t+\t.\tID=id-NC_001133.9:151465..151474;Note=CEN1
NC_001133.9\tRefSeq\tcentromere\t151475\t151557\t.\t+\t.\tID=id-NC_001133.9:151475..151557;Note=CEN1
NC_001133.9\tRefSeq\tcentromere\t151558\t151582\t.\t+\t.\tID=id-NC_001133.9:151558..151582;Note=CEN1
NC_001224.1\tRefSeq\tregion\t1\t85779\t.\t+\t.\tID=NC_001224.1:1..85779;Dbxref=taxon:559292;Name=MT;gbkey=Src;genome=mitochondrion;mol_type=genomic DNA;strain=S288C;Is_circular=true;substrain=FY1679
NC_001224.1\tRefSeq\tgene\t731\t802\t.\t+\t.\tID=gene-tP(UGG)Q;Dbxref=GeneID:854578,SGD:S000007307;Name=tP(UGG)Q;gbkey=Gene;gene_biotype=tRNA
        """
        self.expected_bed_content = """chrI\t640\t700\ttRNA\t+
chrI\t1000\t2000\tgene1(gene)\t+
chrI\t3000\t4000\ttRNA\t-
chrI\t5000\t6000\tTy1(LTR)\t+
"""
        with open(self.gff_file, 'w') as fh:
            fh.write(self.gff_content.strip())
        self.gff_db = create_gff_db(gff_file = self.gff_file)
        assert self.gff_db.count_features_of_type() == 11
        return

    def tearDown(self):
        for test_file in [self.gff_file, self.bed_file]:
            if os.path.exists(test_file):
                os.remove(test_file)
        return super().tearDown()

    def test_parse_gff_to_bed(self):
        from reformat_genome_features_to_bed import parse_gff_to_bed
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
        self.assertEqual(8, self.gff_db.count_features_of_type())
        self.assertEqual(1, self.gff_db.count_features_of_type("centromere"))
        return

    def test_get_buffer_specification_index(self):
        from reformat_genome_features_to_bed import get_buffer_specification_index, BUFFER_SPEC
        test_buffer = BUFFER_SPEC(feature_name="centromere", length=1000, prefix=["peri"], sep="")
        res = get_buffer_specification_index(buffer_spec = [test_buffer])
        for ft_id, features in res.items():
            self.assertEqual(ft_id, "centromere")
            self.assertEqual(1, len(features))

        # # TODO: Add features directly to GFF database
        # self.assertEqual(13, self.gff_db.count_features_of_type())
        return
    
    def test_create_buffered_features(self):
        from reformat_genome_features_to_bed import create_buffered_features
        from gffutils import Feature
        res = create_buffered_features(self.gff_db,
                                       buffer_feature="tRNA",
                                       buffer_length=1000)
        self.assertEqual(['tRNA1', 'tRNA2'], list(res.keys()))
        self.assertEqual(2, len(res['tRNA1']))
        self.assertIsInstance(res['tRNA1'][0], Feature)
        for feature in res['tRNA1']:
            if feature.id == 'tRNA1_upstream':
                self.assertEqual(1, feature.start)
                self.assertEqual(640, feature.end)
            if feature.id == 'tRNA1_downstream':
                self.assertEqual(700, feature.start)
                self.assertEqual(1700, feature.end)

        for feature in res['tRNA2']:
            if feature.id == 'tRNA2_upstream':
                self.assertEqual(4000, feature.start)
                self.assertEqual(5000, feature.end)
            if feature.id == 'tRNA2_downstream':
                self.assertEqual(2000, feature.start)
                self.assertEqual(3000, feature.end)

        # TODO: Add features directly to GFF database
        # self.assertEqual(15, self.gff_db.count_features_of_type())
        # self.assertEqual(2, self.gff_db.count_features_of_type("downstream_tRNA"))
        # self.assertEqual(2, self.gff_db.count_features_of_type("upstream_tRNA"))
        return
    
    def test_load_buffer_feature_specifications(self):
        from reformat_genome_features_to_bed import load_buffer_feature_specifications, BUFFER_SPEC
        buffer_spec_file = "test_buffer_spec.tsv"
        buffer_spec_content = """# Feature name\tBuffer prefix\tSeparator character\tBuffer length (bp)
centromere\tperi\t\t1000
tRNA\tupstream,downstream\t_\t500
"""

        with open(buffer_spec_file, 'w') as fh:
            fh.write(buffer_spec_content.strip())

        buffer_specs = load_buffer_feature_specifications(buffer_spec_file)

        expected_specs = [
            BUFFER_SPEC(feature_name="centromere", prefix=["peri"], sep="", length=1000),
            BUFFER_SPEC(feature_name="tRNA", prefix=["upstream", "downstream"], sep="_", length=500)
        ]

        self.assertEqual(buffer_specs, expected_specs)

        if os.path.exists(buffer_spec_file):
            os.remove(buffer_spec_file)
        return

    def test_get_feature_name(self):
        from reformat_genome_features_to_bed import get_feature_name
        from gffutils import Feature
        from gffutils.feature import feature_from_line
        with self.subTest("Test getting feature name for gene feature"):
            gene_feature = Feature(seqid="chrI", source="Mock", featuretype="gene", start=1000, end=2000, strand="+", id="gene1", attributes={"ID": ["gene1"]})
            self.assertEqual(get_feature_name(gene_feature), "gene1(gene)")

        with self.subTest("Test getting feature name for tRNA feature"):
            trna_feature = Feature(seqid="chrI", source="Mock", featuretype="tRNA", start=640, end=700, strand="+", attributes={"ID": ["tRNA1"]})
            self.assertEqual(get_feature_name(trna_feature), "tRNA")

        with self.subTest("Test getting feature name for mobile_genetic_element feature"):
            repeat_feature = Feature(seqid="NC_001133.9", source="Mock", featuretype="mobile_genetic_element", start=160238, end=166162, strand="-", attributes={"ID": ["id-NC_001133.9:160238..166162"], "gbkey": ["mobile_element"], "mobile_element_type": ["retrotransposon:YARCTy1-1"]})
            self.assertEqual(get_feature_name(repeat_feature), "Ty1(retrotransposon)")

        with self.subTest("Test getting feature name for sequence_feature"):
            etc_feature = feature_from_line("NC_001134.8\tRefSeq\tsequence_feature\t549229\t549251\t.\t+\t.\tID=id-NC_001134.8:549229..549251;Dbxref=SGD:S000077079;Note=ETC8%3B Chromosome-organizing-clamp%3B tethers chromosomal regions to the nuclear periphery%3B binds TFIIIC transcription factor but does not recruit RNA Polymerase III%3B can act as a transcription-blocking insulator or as a heterochromatin barrier element%3B located between divergently transcribed genes RPB5 and CNS1;gbkey=misc_feature")
            self.assertEqual(get_feature_name(etc_feature), "ETC8(misc_feature)")
            nts_feature = feature_from_line("NC_001144.5\tRefSeq\tsequence_feature\t459797\t460711\t.\t-\t.\tID=id-NC_001144.5:459797..460711;Dbxref=SGD:S000029714;Note=NTS1-2%3B Non-transcribed region of the rDNA repeat between the 3' ETS and RDN5%3B required for rDNA repeat expansion;gbkey=misc_feature")
            self.assertEqual(get_feature_name(nts_feature), "NTS1(misc_feature)")

        with self.subTest("Test getting feature name for unknown feature"):
            unknown_feature = Feature(seqid="chrI", source="Mock", featuretype="unknown", start=1000, end=2000, strand="+", attributes={"ID": ["unknown1"]})
            self.assertRaises(ValueError, get_feature_name, unknown_feature)

        with self.subTest("Test parsing Ty1 element information from LTR feature"):
            feature_inst = feature_from_line("NC_001133.9\tRefSeq\tlong_terminal_repeat\t22230\t22552\t.\t+\t.\tID=id-NC_001133.9:22230..22552;Dbxref=SGD:S000006787;Note=Ty1 LTR;gbkey=repeat_region;rpt_type=long_terminal_repeat")
            self.assertEqual(get_feature_name(feature_inst), "Ty1(LTR)")

        with self.subTest("Test getting rRNA feature name"):
            mock_inst = feature_from_line("NC_001144.5\tRefSeq\trRNA\t459676\t459796\t.\t+\t.\tID=rna-NR_132215.1;Parent=gene-YNCL0018W;Dbxref=SGD:S000006479,GeneID:9164945,GenBank:NR_132215.1;Name=NR_132215.1;gbkey=rRNA;gene=RDN5-1;locus_tag=YNCL0018W;orig_transcript_id=gnl|WGS:|mrna.YNCL0018W-RA;product=5S ribosomal RNA;transcript_id=NR_132215.1")
            self.assertEqual(get_feature_name(mock_inst), "5S rRNA")
        return
    

    def test_update_gff_with_additional_annotations(self):
        from reformat_genome_features_to_bed import update_gff_with_additional_annotations, create_gff_db
        annot_table = "test_annotations.bed"
        annot_content = """NC_001133.9\t100\t200\tcustom_feature1\t+
NC_001133.9\t300\t400\tcustom_feature2\t-
"""

        with open(annot_table, 'w') as fh:
            fh.write(annot_content.strip())

        new_features = update_gff_with_additional_annotations(self.gff_db, annot_table)

        self.assertEqual(len(new_features), 4)
        self.assertIn("custom_feature1", new_features)
        self.assertIn("custom_feature2", new_features)
        self.assertIn("subtelomere-NC_001133.9", new_features)

        feature3 = new_features["subtelomere-NC_001133.9"][0]
        self.assertEqual(feature3.chrom, "NC_001133.9")
        self.assertEqual(feature3.start, 1)
        self.assertEqual(feature3.stop, 15000)
        self.assertEqual(feature3.strand, "+")

        feature3 = new_features["subtelomere-NC_001224.1"][1]
        self.assertEqual(feature3.chrom, "NC_001224.1")
        self.assertEqual(feature3.start, 70779)
        self.assertEqual(feature3.stop, 85779)
        self.assertEqual(feature3.strand, "+")

        if os.path.exists(annot_table):
            os.remove(annot_table)
        return
    
    def test_create_subtelomere_features(self):
        from reformat_genome_features_to_bed import create_subtelomere_features
        from gffutils import Feature
        new_features = create_subtelomere_features(self.gff_db, length=1000)
        self.assertEqual(len(new_features), 2)
        self.assertIn("subtelomere-NC_001133.9", new_features)
        self.assertEqual(4, len(new_features["subtelomere-NC_001133.9"]))

        return
    
    def test_main(self):
        from reformat_genome_features_to_bed import main
        main(
            ["--gff_in",
             "../data/output/feature_locations/Saccharomyces_cerevisiae_S288C.gff",
             "--bed_out", 
             "tests/test_output.bed",
             "--buffer_features",
             "../data/input/buffer_regions.tsv"]
        )
        return

if __name__ == "__main__":
    unittest.main()
