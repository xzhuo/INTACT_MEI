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
        self.variant_id = ":".join([self.chrom, str(self.pos), fields[2]])
        self.mei_anno = fields[3]
        if not (self.mei_anno == "none" or self.mei_anno == "noMEI"):
            (self.mei, self.family, self.strand, self.flag) = self.mei_anno.split(':')
        self.haplotypes = set(fields[4].split(','))

    def print(self):
        return self.line

class mergedVariant: # merged variants at the same coordinates.
    def __init__(self, variant):
        self.chrom = variant.chrom
        self.variant_ids = []
        self.mei = {}
        self.add_variant(variant)

    def add_variant(self, variant):
        self.pos = variant.pos
        self.variant_ids.append(variant.variant_id)
        if variant.mei_anno == "none":
            if hasattr(self, 'missing_haps'):
                self.missing_haps = self.missing_haps.union(variant.haplotypes)
            else:
                self.missing_haps = variant.haplotypes
        elif variant.mei_anno == "noMEI":
            if hasattr(self, 'empty_haps'):
                self.empty_haps = self.empty_haps.intersection(variant.haplotypes)
            else:
                self.empty_haps = variant.haplotypes
        else:
            if self.mei.get(variant.family):
                if self.mei[variant.family].get(variant.strand):
                    self.mei[variant.family][variant.strand]['meis'].append(variant.mei)
                    self.mei[variant.family][variant.strand]['intact'].append(variant.flag)
                    self.mei[variant.family][variant.strand]['haps'] = self.mei[variant.family][variant.strand]['haps'].union(variant.haplotypes)
                else:
                    self.mei[variant.family][variant.strand] = {'meis':[variant.mei], 'intact':[variant.flag], 'haps':variant.haplotypes}
            else:
                self.mei[variant.family] = {variant.strand:{'meis':[variant.mei], 'intact':[variant.flag], 'haps':variant.haplotypes}}

    def combine_variants(self,merged_variant):
        pass

    def print(self):
        samples_str = ','.join(sorted(self.haplotypes))
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.variant_ids}\t{samples_str}"

def merge_per_pos(input_tsv):
    list_merged = []
    with open(input_tsv, 'r') as f:
        for line in f:
            variant = Variant(line)
            if list_merged and list_merged[-1].chrom == variant.chrom and list_merged[-1].end == variant.pos:
                list_merged[-1].add_variant(variant)
            else:
                list_merged.append(mergedVariant(variant))
    return list_merged


def main():
    p = argparse.ArgumentParser(description="merge MEI calling variants from MCgraph VCF file.")
    p.add_argument("--input", required=True, help="Input TSV file")
    p.add_argument("--out", dest="out_tsv", required=True, help="Output intact-only TSV path")
    args = p.parse_args()
    input_tsv = args.input
    out_tsv = args.out_tsv
    list_merged = merge_per_pos(input_tsv)
    for variant in list_merged:
        chrom = variant.chrom
        pos = variant.pos
        for family in variant.mei:
            for strand in variant.mei[family]:
                details = variant.mei[family][strand]
                all_dict[family][strand] = details

if __name__ == "__main__":
    main()