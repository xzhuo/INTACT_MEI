VCF = "HG005.pav.insdel.vcf.gz"
HG38 = "hg38.fa"

rule all:
  input:
    "HG005.indel.fa",
    "HG005.indel.fa.out",
    "HG005.limeaid.v134.tsv",
    "HG005.limeaid.v134.intact.tsv",
    "HG005.info.hg38.txt.gz",
    "HG005.pav.intact_mei.hg38.vcf.gz"

rule annotate_info:
    # bcftools not in the docker image.
    input:
        vcf = VCF,
        info = "HG005.info.hg38.txt.gz",
        header = "scripts/intact.info.header.txt"
    output:
        "HG005.pav.intact_mei.hg38.vcf.gz"
    shell:
        "bcftools annotate -a {input.info} -c CHROM,POS,~ID,INFO/INTACT_MEI,INFO/INTACT_FLAG -h {input.header} -Oz -o {output} {input.vcf}"

rule zip_info:
    # bgzip and tabix not in the docker image.
    input:
        "HG005.info.hg38.txt"
    output:
        "HG005.info.hg38.txt.gz"
    shell:
        "bgzip {input}; tabix -s 1 -b 2 -e 2 {input}.gz"

rule intact_mei:
    input:
        limeaid = "HG005.limeaid.v134.tsv",
        vcf = VCF
    output:
        tsv = "HG005.limeaid.v134.intact.tsv",
        info = "HG005.info.hg38.txt"
    shell:
        "python3 scripts/limeaid_intact_mei.py --input {input.limeaid} --vcf {input.vcf} --flag-table {output.info} --out {output.tsv}"

rule limeaid:
    input:
        fa = "HG005.indel.fa",
        rmsk = "HG005.indel.fa.out"
    output:
        "HG005.limeaid.v134.tsv"
    shell:
        "python3 /opt/src/L1ME-AID/limeaid.v1.3.4-beta.py -i {input.fa} -g {HG38} -r {input.rmsk} -o {output}"

rule rmsk:
    # mount famdb to /opt/RepeatMasker/Libraries/famdb
    input:
        "HG005.indel.fa"
    output:
        "HG005.indel.fa.out"
    shell:
        "/opt/RepeatMasker/RepeatMasker -species human {input}"

rule extract_fa:
    input:
        VCF
    output:
        "HG005.indel.fa"
    shell:
        "perl scripts/extract_fa.pl {input} > {output}"