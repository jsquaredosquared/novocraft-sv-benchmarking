configfile: "../../config/config.yaml"

SAMPLES = config["samples"]
TRUTH_SETS = config["truth_sets"]
REFERENCE = config["reference"]
ALIGNERS = config["aligners"]
CALLERS = config["callers"]

CATEGORIES=["ALL", "DEL", "INS"]


rule all:
    input:
        expand("../../outputs/truvari/{sample}.{aligner}.{caller}.{category}.truvari-bench.json", sample=SAMPLES, aligner=ALIGNERS, caller=CALLERS, category=CATEGORIES)


def get_filter_expression(wildcards):
    match wildcards.category.split('-'):
        case [category]:
            return 'True' if category == "ALL" else f'INFO["SVTYPE"] == "{category}"'
        case [min_len, max_len]:
            return f'{int(min_len)}<=INFO["SVLEN"]<={int(max_len)}'


rule split_vcf_into_category:
    input:
        "../../outputs/{caller}/{sample}.{aligner}.{caller}.vcf.gz"
    output:
        vcf_file = temp("../../outputs/{caller}/{sample}.{aligner}.{caller}.{category}.vcf.gz"),
        index_file = temp("../../outputs/{caller}/{sample}.{aligner}.{caller}.{category}.vcf.gz.tbi")
    conda:
        "../envs/vembrane.yaml"
    log:
        "../../logs/{sample}.{aligner}.{caller}.split_vcf_{category}.log"
    params:
        expression = get_filter_expression,
        options = lambda wildcards: "--overwrite-number-info SVLEN=1" if wildcards.caller=="manta" else ""
    shell:
        "(vembrane filter {params.options} '{params.expression}' {input} "
        "| bcftools sort - "
        "| bgzip --stdout - "
        "> {output.vcf_file} "
        "&& tabix -p vcf {output.vcf_file} "
        ")2> {log}"


rule compare_to_truth_set:
    input:
        multiext("../../outputs/{caller}/{sample}.{aligner}.{caller}.{category}.vcf", ".gz", ".gz.tbi"),
    output:
        "../../outputs/truvari/{sample}.{aligner}.{caller}.{category}.truvari-bench.json",
    conda:
        "../envs/truvari.yaml",
    log:
        "../../logs/{sample}.{aligner}.{caller}.{category}.truvari-bench.log",
    params:
        truth_set=lambda wildcards: TRUTH_SETS[wildcards.sample].replace(".vcf", f".{wildcards.category}.vcf"),
        out_dir=lambda wildcards: f"../../outputs/truvari/{wildcards.sample}.{wildcards.aligner}.{wildcards.caller}.{wildcards.category}/",
        regions="../../resources/sv-benchmarks/HG002/HG002_SVs_Tier1_v0.6.bed",
    shell:
        "(truvari bench "
        "--base {params.truth_set} "
        "--comp {input[0]} "
        "--output {params.out_dir} "
        "--includebed {params.regions} "
        "--pctseq 0 "
        "--sizemin 50 "
        "--sizemax 1_000_000 "
        "--passonly "
        "&& cat {params.out_dir}/summary.json > {output} "
        ")2> {log}"
