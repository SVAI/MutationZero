#!/usr/bin/env python
"""
Takes a tsv file (CADD output) and a vcf from 1KGP (1000 genomes project)
and filters out the 1KGP variants from the variants scored with CADD.
"""

import vcf
import sys

infile = sys.argv[1]
outfile = sys.argv[2]
filtered_cadd = sys.argv[3]

reader_1000g = vcf.Reader(open(infile, 'r'))

common_indels = {} 

# gather all 1000g indels
i = 0
for record in reader_1000g:
    i += 1
    chr_pos = (str(record.CHROM), int(record.POS))
    common_indels[chr_pos] = (record.REF, record.ALT[0])
    if i % 10000 == 0:
        print i

writer = open(outfile, 'w')
removed = 0
# check if there are any identical mutations in the same position
with open(filtered_cadd, 'r') as cadd:
    next(cadd)
    next(cadd)
    for line in cadd:
        spl = line.split()
        chr_pos = (str(spl[0]), int(spl[1]))
        if chr_pos in common_indels and common_indels[chr_pos][0] == spl[2] and common_indels[chr_pos][1]  == spl[3]:
            removed += 1
            continue
        else:
            writer.write(line)
writer.close()

print len("We removed " + str(removed))
