FILES = glob_wildcards('{s}.vcf.gz')
VCFs = FILES.s
HG38 = "hg38.fa"

rule all:
  input:
    expand("{s}.indel.fa", s=VCFs),
    expand("{s}.indel.fa.out", s=VCFs),
    expand("{s}.limeaid.v134.tsv", s=VCFs),
    expand("{s}.limeaid.v134.intact.tsv", s=VCFs),
    # expand("{s}.info.hg38.txt.gz", s=VCFs),
    # expand("{s}.info.hg38.txt.gz.tbi", s=VCFs),
    expand("{s}.intact_mei.hg38.vcf.gz", s=VCFs),
    expand("{s}.intact_mei.hg38.vcf.gz.tbi", s=VCFs),

rule annotate_info:
    # bcftools not in the dockerfile.
    input:
        vcf = "{s}.vcf.gz",
        info = "{s}.info.hg38.txt.gz",
        header = "scripts/intact.info.header.txt"
    output:
        gz = "{s}.intact_mei.hg38.vcf.gz",
        tbi = "{s}.intact_mei.hg38.vcf.gz.tbi"
    container:
        "docker://xiaoyuz/biotools:latest"
    shell:
        "bcftools annotate -a {input.info} -c CHROM,POS,~ID,INFO/INTACT_MEI,INFO/INTACT_FLAG -h {input.header} -Oz -o {output.gz} {input.vcf} && bcftools index -t {output.gz}"

rule zip_info:
    # bgzip and tabix not in the dockerfile.
    input:
        "{s}.info.hg38.txt"
    output:
        gz = "{s}.info.hg38.txt.gz",
        tbi = "{s}.info.hg38.txt.gz.tbi"
    container:
        "docker://xiaoyuz/biotools:latest"
    shell:
        "bgzip {input} && tabix -s 1 -b 2 -e 2 {input}.gz"

rule intact_mei:
    input:
        limeaid = "{s}.limeaid.v134.tsv",
        vcf = "{s}.vcf.gz"
    output:
        tsv = "{s}.limeaid.v134.intact.tsv",
        info = "{s}.info.hg38.txt"
    container:
        "docker://xiaoyuz/biotools:latest"
    shell:
        "python3 scripts/limeaid_intact_mei.py --input {input.limeaid} --vcf {input.vcf} --flag-table {output.info} --out {output.tsv}"

rule limeaid:
    input:
        fa = "{s}.indel.fa",
        rmsk = "{s}.indel.fa.out"
    output:
        "{s}.limeaid.v134.tsv"
    container:
        "docker://xiaoyuz/l1me-aid:1.3.4"
    shell:
        "python3 /opt/src/L1ME-AID/limeaid.v1.3.4-beta.py -i {input.fa} -g {HG38} -r {input.rmsk} -o {output}"
        # "python3 /opt/src/L1ME-AID/limeaid.v1.3.4-beta.py -i {input.fa} -r {input.rmsk} -o {output}" # much faster without TSD

rule rmsk:
    # mount famdb to /opt/RepeatMasker/Libraries/famdb
    input:
        "{s}.indel.fa"
    output:
        "{s}.indel.fa.out"
    container:
        "docker://xiaoyuz/l1me-aid:1.3.4"
    shell:
        "/opt/RepeatMasker/RepeatMasker -species human {input}"

rule extract_fa:
    input:
        "{s}.vcf.gz"
    output:
        "{s}.indel.fa"
    container:
        "docker://xiaoyuz/biotools:latest"
    shell:
        "perl scripts/extract.indel.seq.pl {input} {output}"