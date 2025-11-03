# merge MEI calling variants from MCgraph VCF file
# input file lines: 
'''
chr10   296371  .       none    NA20762_hap2
chr10   296371  1       AluY:SINE/Alu:+:INTACT  NA19159_hap1
chr10   331885  0       noMEI   HG02922_hap1,HG02922_hap2
'''

import sys
import argparse

class Variant:
    def __init__(self, line):
        self.line = line.strip()
        self.process_line()

    def process_line(self):
        fields = self.line.split()
        self.chrom = fields[0]
        self.pos = int(fields[1])
        self.variant_id = fields[2] # don't int it, could be '.'
        self.mei_anno = fields[3]
        if not (self.mei_anno == "none" or self.mei_anno == "noMEI"):
            (self.mei, self.family, self.strand, self.flag) = self.mei_anno.split(':')
        self.samples = fields[4].split(',')

    def print(self):
        return self.line

class mergedVariant:
    def __init__(self, variant):
        self.chrom = variant.chrom
        self.start = variant.pos
        self.end = variant.pos
        match variant.mei_anno:
            case "none":
                self.variant_num = 1
            case "noMEI":
                self.variant_num = 2
            case _:
                self.variant_num = 0

    def add_variant(self, variant):
        self.flags.add(variant.flag)
        for sample in variant.samples:
            self.samples.add(sample)

    def print(self):
        flags_str = ','.join(sorted(self.flags))
        samples_str = ','.join(sorted(self.samples))
        return f"{self.chrom}\t{self.pos}\t{self.variant_num}\t{self.mei}:{self.family}:{self.strand}:{flags_str}\t{samples_str}"

def main():
    p = argparse.ArgumentParser(description="merge MEI calling variants from MCgraph VCF file.")
    p.add_argument("--input", required=True, help="Input TSV file")
    p.add_argument("--out", dest="out_tsv", required=True, help="Output intact-only TSV path")
    args = p.parse_args()
    input_tsv = args.input
    out_tsv = args.out_tsv
    pass

if __name__ == "__main__":
    main()