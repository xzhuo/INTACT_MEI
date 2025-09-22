FILES = glob_wildcards('input/{s}.vcf.gz') # leave input files in the input folder. {s} will be the sample name prefix.
VCFs = FILES.s
HG38 = "hg38.fa"

rule all:
  input:
    # direct output to the output folder
    expand("output/{s}.indel.fa", s=VCFs), 
    expand("output/{s}.indel.fa.out", s=VCFs),
    expand("output/{s}.limeaid.v134.tsv", s=VCFs),
    expand("output/{s}.limeaid.v134.intact.tsv", s=VCFs),
    # expand("output/{s}.info.hg38.txt.gz", s=VCFs),
    # expand("output/{s}.info.hg38.txt.gz.tbi", s=VCFs),
    expand("output/{s}.intact_mei.hg38.vcf.gz", s=VCFs),
    expand("output/{s}.intact_mei.hg38.vcf.gz.tbi", s=VCFs)

rule annotate_info:
    # bcftools not in the dockerfile.
    input:
        vcf = "input/{s}.vcf.gz",
        info = "output/{s}.info.hg38.txt.gz",
        header = "scripts/intact.info.header.txt"
    output:
        gz = "output/{s}.intact_mei.hg38.vcf.gz",
        tbi = "output/{s}.intact_mei.hg38.vcf.gz.tbi"
    container:
        "docker://xiaoyuz/biotools:latest"
    shell:
        "bcftools annotate -a {input.info} -c CHROM,POS,~ID,INFO/INTACT_MEI,INFO/INTACT_FLAG -h {input.header} -Oz -o {output.gz} {input.vcf} && bcftools index -t {output.gz}"

rule zip_info:
    # bgzip and tabix not in the dockerfile.
    input:
        "output/{s}.info.hg38.txt"
    output:
        gz = "output/{s}.info.hg38.txt.gz",
        tbi = "output/{s}.info.hg38.txt.gz.tbi"
    container:
        "docker://xiaoyuz/biotools:latest"
    shell:
        "bgzip {input} && tabix -s 1 -b 2 -e 2 {input}.gz"

rule intact_mei:
    input:
        limeaid = "output/{s}.limeaid.v134.tsv",
        vcf = "input/{s}.vcf.gz"
    output:
        tsv = "output/{s}.limeaid.v134.intact.tsv",
        info = "output/{s}.info.hg38.txt"
    container:
        "docker://xiaoyuz/biotools:latest"
    shell:
        "python3 scripts/limeaid_intact_mei.py --input {input.limeaid} --vcf {input.vcf} --flag-table {output.info} --out {output.tsv}"

rule limeaid:
    input:
        fa = "output/{s}.indel.fa",
        rmsk = "output/{s}.indel.fa.out"
    output:
        "output/{s}.limeaid.v134.tsv"
    container:
        "docker://xiaoyuz/l1me-aid:1.3.4"
    shell:
        "python3 /opt/src/L1ME-AID/limeaid.v1.3.4-beta.py -i {input.fa} -g {HG38} -r {input.rmsk} -o {output}"
        # "python3 /opt/src/L1ME-AID/limeaid.v1.3.4-beta.py -i {input.fa} -r {input.rmsk} -o {output}" # much faster without TSD

rule rmsk:
    # mount famdb to /opt/RepeatMasker/Libraries/famdb
    input:
        "output/{s}.indel.fa"
    output:
        "output/{s}.indel.fa.out"
    container:
        "docker://xiaoyuz/l1me-aid:1.3.4"
    shell:
        "/opt/RepeatMasker/RepeatMasker -species human {input}"

rule extract_fa:
    input:
        "input/{s}.vcf.gz"
    output:
        "output/{s}.indel.fa"
    container:
        "docker://xiaoyuz/biotools:latest"
    shell:
        "perl scripts/extract.indel.seq.pl {input} {output}"