# LIMEaid Intact MEI Pipeline

End-to-end workflow to filter LIMEaid TSVs for intact/LTR events and annotate a VCF with INTACT flags using bcftools.

## Components

- `limeaid_filter.py` — Produce:
- Intact/LTR-only TSV (`--out`)
- 5-column flag table (`--flag-table`): `CHROM POS ID INTACT_MEI INTACT_FLAG`
- CHROM/POS are read from the given VCF by matching variant `ID` (no ID parsing).
- `annotate_intact_info.sh` — Adds `INTACT_MEI` and `INTACT_FLAG` to VCF INFO using the flag table (CHROM,POS,~ID match).
- `limeaid_process.sh` — One-shot script: filtering + annotation.
- `environment.yml` — Conda environment for reproducibility.

## Requirements

- Python 3.10+
- bcftools, bgzip, tabix (htslib)

## Quick Start

1) Generate outputs from a LIMEaid TSV and VCF:
   - `./limeaid_filter.py --input input.limeaid.tsv --vcf input.vcf.gz --flag-table output.info.txt --out output.intact.tsv`
2) Annotate the VCF
   - One-shot:
     - `./limeaid_process.sh input.limeaid.tsv input.vcf.gz output.intact.tsv output.info.txt output.annot.vcf.gz`
   - Manual:
     - `sort -k1,1 -k2,2n output.info.txt | bgzip -c > output.info.txt.gz`
     - `tabix -s 1 -b 2 -e 2 output.info.txt.gz`
     - `printf '%s\n' '##INFO=<ID=INTACT_MEI,Number=1,Type=String,Description="MEI subtype from limeaid intact flags">' '##INFO=<ID=INTACT_FLAG,Number=1,Type=String,Description="INTACT status from limeaid (e.g., INTACT, INTACT_3end)">' > intact.info.hdr`
     - `bcftools annotate -a output.info.txt.gz -c CHROM,POS,~ID,INFO/INTACT_MEI,INFO/INTACT_FLAG -h intact.info.hdr -Oz -o output.annot.vcf.gz input.vcf.gz`
     - `bcftools index -t output.annot.vcf.gz`

## Behavior

- Skips header lines and rows with empty class (column 9) in the LIMEaid TSV.
- Applies the same flagging as the Perl script (writes to column 11): `soloLTR`, `ERV-LTR`, `LTR-ERV-LTR`, `INTACT`, `INTACT_3end`.
- Flag table has no header; CHROM/POS come from the VCF via ID match. Rows without a matching ID in the VCF are omitted.

## Review

- Header lines: `bcftools view -h output.annot.vcf.gz | grep -E 'INTACT_(MEI|FLAG)'`
- Annotated rows: `bcftools view -H -i 'INFO/INTACT_FLAG!=""' output.annot.vcf.gz | head`
- Counts: `bcftools view -H -i 'INFO/INTACT_FLAG!=""' output.annot.vcf.gz | wc -l`

