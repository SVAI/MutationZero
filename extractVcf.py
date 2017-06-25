#!/usr/bin/env python
import vcf
import sys

"""
Takes as input a VEP annotated vcf file and extracts
only variants associated with (our) genes of interest
"""
infile = sys.argv[1]
outfile = sys.argv[2]


class RangeDict():
    def __init__(self):
        self._dict = {}

    def __getitem__(self, key):
        for k, v in self._dict.items():
            if k[0] <= key < k[1]:
                return v
        raise KeyError("Key not found!")

    def __setitem__(self, key, value):
        if len(key) == 2:
            if key[0] < key[1]:
                self._dict.__setitem__((key[0], key[1]), value)

    def __contains__(self, key):
        try:
            return bool(self.__getitem__(key))
        except KeyError:
            return False


reader = vcf.Reader(open(infile, 'r'))
vcf_writer = vcf.Writer(open(outfile, 'w'), reader)

def flatten(li):
    return [item for sublist in li for item in sublist]

def get_genes(gene_info):
    """
    Returns a list of tuples (gene_symbol, gene_id)
    """
    return [i.split(":") for i in gene_info]

def write_records(vcf_writer, list_of_records):
    for r in list_of_records:
        vcf_writer.write_record(r)

# EGFR: chr7:55247443-55256642
# RPS3A: chr4:152020725-152025804
# RPS3: chr11:75110535-75117957
# RPS2: chr16:2012062-2014827
# RPS7: chr2:3622853-3628509
# RPS5: chr19:58898636-58906171
# RPS4X: chrX:71492453-71497141
# RPL14: chr3:40498783-40503863
# RPL4: chr15:66791653-66797193
# RPL18A: chr19:17970687-17974133
# HDAC1: chr1:32757708-32799224
# RPL18: chr19:49118584-49122675
chroms = ['7:55247443-55256642', '4:152020725-152025804', '11:75110535-75117957', '16:2012062-2014827', '2:3622853-3628509',
          '19:58898636-58906171', 'X:71492453-71497141', '3:40498783-40503863', '15:66791653-66797193', '19:17970687-17974133', '1:32757708-32799224',
          '19:49118584-49122675', '22:29999545-30094589']
r_genes = RangeDict()
for i in chroms:
    first = i.split(':')[1].split('-')[0]
    second = i.split(':')[1].split('-')[1]
    r_genes[(int(first),str(second))] = i.split(':')[0]

l = set()
last_header = "Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|EXON|INTRON|HGNC|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|DISTANCE|CLIN_SIG|CANONICAL|SIFT|PolyPhen|GMAF|BIOTYPE|ENSP|DOMAINS|CCDS|HGVSc|HGVSp|CELL_TYPE".split("|")
for record in reader:
    # Using chrom and pos to filter by
    chrom = str(record.CHROM)
    pos = int(record.POS)
    info = record.INFO

    if 'CSQ' in info:
        data = info['CSQ'][0].split('|')
        for indx, col in enumerate(last_header):
            record.add_info("added."+col, data[indx])

    if pos in r_genes and r_genes[pos] == chrom: #filter by chrom and pos
        l.add(record)

write_records(vcf_writer, l)
