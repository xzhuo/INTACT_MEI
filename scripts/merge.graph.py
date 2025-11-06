# merge MEI calling variants from MCgraph VCF file
# input file lines: 
'''
chr10   296371  .       none    NA20762_hap2
chr10   296371  1       AluY:SINE/Alu:+:INTACT  NA19159_hap1
chr10   331885  0       noMEI   HG02922_hap1,HG02922_hap2
'''

import sys
import argparse
import statistics
from statistics import mode

class Variant:
    def __init__(self, line):
        self.line = line.strip()
        self.process_line()

    def process_line(self):
        fields = self.line.split()
        self.chrom = fields[0]
        self.pos = int(fields[1])
        self.is_ref = True if fields[2] == '0' else False
        self.variant_id = ":".join([self.chrom, str(self.pos), fields[2]])
        self.mei_anno = fields[3]
        if not (self.mei_anno == "none" or self.mei_anno == "noMEI"):
            (self.mei, self.family, self.strand, self.flag) = self.mei_anno.split(':')
        self.length = int(fields[4])
        self.haplotypes = set(fields[5].split(','))

    def print(self):
        return self.line

class mergedVariant: # merged variants at the same coordinates.
    def __init__(self, variant):
        self.chrom = variant.chrom
        self.start = variant.pos
        self.variant_ids = []
        self.missing_haps = set()
        self.empty_haps = set()
        self.mei = {}
        self.add_variant(variant)

    def indel_validate(self):
        if hasattr(self, 'ins_haps') and len(self.ins_haps) > 0:
            if len(self.mei) <= 1:
                for family in self.mei:
                    for strand in self.mei[family]:
                        self.mei[family][strand]['haps'].update(self.ins_haps)
                del self.ins_haps
            else:
                breakpoint()
        if not hasattr(self, 'indel') or not self.indel:
            if (len(self.empty_haps) == 0):
                self.indel = "INS"
            elif len(self.mei) == 0:
                self.indel = "DEL"
            else:
                self.indel = "UNKNOWN"

    def add_variant(self, variant):
        self.end = variant.pos
        self.variant_ids.append(variant.variant_id)
        if variant.is_ref:
            if variant.mei_anno == "noMEI" and variant.length == 0:
                self.indel = "INS"
            else:
                self.indel = "DEL"
        if variant.mei_anno == "none":
            if len(self.missing_haps) > 0:
                self.missing_haps = self.missing_haps.union(variant.haplotypes)
            else:
                self.missing_haps = variant.haplotypes
        elif variant.mei_anno == "noMEI":
            if variant.length == 0:
                if len(self.empty_haps) > 0:
                    self.empty_haps = self.empty_haps.intersection(variant.haplotypes)
                else:
                    self.empty_haps = variant.haplotypes
            else:
                if not hasattr(self, 'ins_haps'):
                    self.ins_haps = set()
                self.ins_haps = self.ins_haps.union(variant.haplotypes)
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
        self.start = min(self.start, merged_variant.start)
        self.end = max(self.end, merged_variant.end)
        self.variant_ids.extend(merged_variant.variant_ids)
        self.missing_haps = self.missing_haps.union(merged_variant.missing_haps)
        self.empty_haps = self.empty_haps.intersection(merged_variant.empty_haps)
        if hasattr(merged_variant, 'indel') and hasattr(self, 'indel'):
            if self.indel != merged_variant.indel:
                self.indel = 'UNKNOWN'
        elif hasattr(merged_variant, 'indel'):
            self.indel = merged_variant.indel
        for family in merged_variant.mei:
            for strand in merged_variant.mei[family]:
                details = merged_variant.mei[family][strand]
                if self.mei.get(family):
                    if self.mei[family].get(strand):
                        self.mei[family][strand]['meis'].extend(details['meis'])
                        self.mei[family][strand]['intact'].extend(details['intact'])
                        self.mei[family][strand]['haps'] = self.mei[family][strand]['haps'].union(details['haps'])
                    else:
                        self.mei[family][strand] = details
                else:
                    self.mei[family] = {strand:details}

    def print(self):
        samples_list = []
        if len(self.empty_haps) > 0:
            empty_haps_str = ','.join(self.empty_haps)
            samples_list.append(f"noMEI={empty_haps_str}")
        for family in self.mei:
            for strand in self.mei[family]:
                details = self.mei[family][strand]
                the_mei = mode(details['meis'])
                the_intact = mode(details['intact'])
                key = f"{the_mei}:{family}:{strand}:{the_intact}"
                value = ','.join(details['haps'])
                mei_str = f"{key}={value}"
                samples_list.append(mei_str)
        samples_str = ';'.join(samples_list)
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.indel}\t{self.variant_ids}\t{samples_str}"

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
    list_combied_merged = []
    for merged_variant in list_merged:
        merged_variant.indel_validate()

        chrom = merged_variant.chrom
        pos = merged_variant.start
        if list_combied_merged and list_combied_merged[-1].chrom == chrom and list_combied_merged[-1].end + 100 >= pos:
            to_be_combined = False
            for family in merged_variant.mei:
                for strand in merged_variant.mei[family]:
                    if list_combied_merged[-1].mei.get(family):
                        if list_combied_merged[-1].mei[family].get(strand):
                            to_be_combined = True
            if to_be_combined:
                list_combied_merged[-1].combine_variants(merged_variant)
            else:
                list_combied_merged.append(merged_variant)
        else:
            list_combied_merged.append(merged_variant)

    with open(out_tsv, 'w') as f:
        for merged_variant in list_combied_merged:
            f.write(merged_variant.print() + '\n')

if __name__ == "__main__":
    main()