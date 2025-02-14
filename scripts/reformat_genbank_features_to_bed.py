#!/usr/bin/env python
# encoding: utf-8
"""
bed_from_genbank.py

grab the gene records from a genbank file (edit for other record types).

- requires:  biopython

"""

from Bio import SeqIO
import sys
import pdb

def main():
    outf = open(sys.argv[2], 'w')
    for record in SeqIO.parse(open(sys.argv[1], "rU"), "genbank") :
        for feature in record.features:
            if feature.type == 'gene':
                start = feature.location.start.position
                stop = feature.location.end.position
                try:
                    #name = "gene"
                    name = feature.qualifiers['gene'][0] + "(gene)"
                except:
                    # some features only have a locus tag
                    name = feature.qualifiers['locus_tag'][0] + "(gene)"
                if feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
                outf.write(bed_line)
            if feature.type == 'telomere':
                start = feature.location.start.position
                stop = feature.location.end.position
                try:
                    name = "telomere"
                except:
                    name = "false"
                if feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
                outf.write(bed_line)
            if feature.type == 'centromere':
                start = feature.location.start.position
                stop = feature.location.end.position
                try:
                    name = "centromere"
                except:
                    name = "false"
                if feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
                outf.write(bed_line)
            if feature.type == 'rep_origin':
                start = feature.location.start.position
                stop = feature.location.end.position
                try:
                    name = "rep_origin"
                except:
                    name = "false"
                if feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
                outf.write(bed_line)
            if feature.type == 'tRNA':
                start = feature.location.start.position
                stop = feature.location.end.position
                try:
                    name = "tRNA"
                except:
                    name = "false"
                if feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
                outf.write(bed_line)
            if feature.type == 'repeat_region':
                start = feature.location.start.position
                stop = feature.location.end.position
                try:
                    name = "LTR"
                except:
                    name = "false"
                if feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
                outf.write(bed_line)
            if feature.type == 'mobile_element':
                start = feature.location.start.position
                stop = feature.location.end.position
                try:
                    name = "LTR"
                except:
                    name = "false"
                if feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                bed_line = "{0}\t{1}\t{2}\t{3}\n".format(start, stop, name, strand)
                outf.write(bed_line)
        outf.close()

if __name__ == '__main__':
    main()
