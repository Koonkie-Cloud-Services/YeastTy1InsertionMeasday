#!/usr/bin/env python

import sys
import re

control = [
    (1000, 'M'),
    (900, 'CM'),
    (500, 'D'),
    (400, 'CD'),
    (100, 'C'),
    (90, 'XC'),
    (50, 'L'),
    (40, 'XL'),
    (10, 'X'),
    (9, 'IX'),
    (5, 'V'),
    (4, 'IV'),
    (1, 'I')
]


def int_to_roman(num):
    roman_num = []
    for x, c in control:
        if (r := num // x):
            roman_num.append(f'{c*r}')
            if (num := num - r * x) == 0:
                break
    return ''.join(roman_num)

def main():
    fh_out = open(sys.argv[2], 'w')

    with open(sys.argv[1]) as fh:
        for line in fh:
            row = line.strip().split("\t")
            # Skip head if present
            if row[0] == "chr":
                pass
            else:

                try:
                    chrom_num = int(re.sub("chr", "", row[0]))
                    row[0] = "chr" + int_to_roman(chrom_num)
                except ValueError:
                    pass
            fh_out.write("\t".join(row) + "\n")

    fh_out.close()
    return

if __name__ == "__main__":
    main()
