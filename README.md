# LIMEaid Intact MEI Pipeline

End-to-end workflow to filter LIMEaid TSVs for intact/LTR events and annotate a VCF with INTACT flags using bcftools.

## Components

- `limeaid_filter.py` — Produce:
  - Intact/LTR-only TSV (`--out`)
  - 5-column flag table (`--flag-table`): `CHROM POS ID INTACT_MEI INTACT_FLAG`
  - CHROM/POS are read from the given VCF by matching variant `ID` (no ID parsing).

## Requirements

- Python 3.10+
- bcftools, bgzip, tabix (htslib)

## Quick Start

1) Generate outputs from a LIMEaid TSV and VCF:
   - `./limeaid_filter.py --input input.limeaid.tsv --vcf input.vcf.gz --flag-table output.info.txt --out output.intact.tsv`
2) Annotate the VCF
   - Manual:
     - `sort -k1,1 -k2,2n output.info.txt | bgzip -c > output.info.txt.gz`
     - `tabix -s 1 -b 2 -e 2 output.info.txt.gz`
     - `bcftools annotate -a output.info.txt.gz -c CHROM,POS,~ID,INFO/INTACT_MEI,INFO/INTACT_FLAG -h scripts/intact.info.header.txt -Oz -o output.annot.vcf.gz input.vcf.gz`
     - `bcftools index -t output.annot.vcf.gz`

## Behavior

- Skips header lines and rows with empty class (column 9) in the LIMEaid TSV.
- Applies the same flagging as the Perl script (writes to column 11): `soloLTR`, `ERV-LTR`, `LTR-ERV-LTR`, `INTACT`, `INTACT_3end`.
- Flag table has no header; CHROM/POS come from the VCF via ID match. Rows without a matching ID in the VCF are omitted.

## Review

- Header lines: `bcftools view -h output.annot.vcf.gz | grep -E 'INTACT_(MEI|FLAG)'`
- Annotated rows: `bcftools view -H -i 'INFO/INTACT_FLAG!=""' output.annot.vcf.gz | head`
- Counts: `bcftools view -H -i 'INFO/INTACT_FLAG!=""' output.annot.vcf.gz | wc -l`
