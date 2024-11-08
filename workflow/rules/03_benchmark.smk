def get_filter_expression(wildcards):
    match wildcards.svtype.split("-"):
        case [svtype]:
            return "True" if svtype == "ALL" else f'INFO["SVTYPE"] == "{svtype}"'
        case [min_len, max_len]:
            return f'{int(min_len)}<=INFO["SVLEN"]<={int(max_len)}'


rule split_truth_set_by_svtype:
    input:
        truth_set=lambda wildcards: config["truth_sets"][wildcards.sample],
    output:
        multiext(
            "resources/sv-benchmarks/{sample}/{sample}.{svtype}_truth-set.vcf",
            ".gz",
            ".gz.tbi",
        ),
    conda:
        "../envs/vembrane.yaml"
    params:
        expression=get_filter_expression,
    log:
        "logs/{sample}.split_truth_set_{svtype}.log"
    shell:
        "(vembrane filter '{params.expression}' {input} "
        "| bgzip --stdout - "
        "> {output[0]} "
        "&& tabix -p vcf {output[0]} "
        ")2> {log}"
        

rule split_vcf_by_svtype:
    input:
        "outputs/{caller}/{sample}.{aligner}.{caller}.vcf.gz",
    output:
        vcf_file=temp("outputs/{caller}/{sample}.{aligner}.{caller}.{svtype}.vcf.gz"),
        index_file=temp(
            "outputs/{caller}/{sample}.{aligner}.{caller}.{svtype}.vcf.gz.tbi"
        ),
    conda:
        "../envs/vembrane.yaml"
    log:
        "logs/{sample}.{aligner}.{caller}.split_vcf_{svtype}.log",
    params:
        expression=get_filter_expression,
        options=lambda wildcards: (
            "--overwrite-number-info SVLEN=1" if wildcards.caller == "manta" else ""
        ),
    shell:
        "(vembrane filter {params.options} '{params.expression}' {input} "
        "| bcftools sort - "
        "| bgzip --stdout - "
        "> {output.vcf_file} "
        "&& tabix -p vcf {output.vcf_file} "
        ")2> {log}"


rule compare_to_truth_set:
    input:
        comp_file="outputs/{caller}/{sample}.{aligner}.{caller}.{svtype}.vcf.gz",
        index_file="outputs/{caller}/{sample}.{aligner}.{caller}.{svtype}.vcf.gz.tbi",
        truth_set="resources/sv-benchmarks/{sample}/{sample}.{svtype}_truth-set.vcf.gz",
        truth_set_regions=lambda wildcards: config["truth_set_regions"][wildcards.sample]
    output:
        "outputs/truvari/{sample}.{aligner}.{caller}.{svtype}.truvari-bench.json",
    conda:
        "../envs/truvari.yaml"
    log:
        "logs/{sample}.{aligner}.{caller}.{svtype}.truvari-bench.log",
    params:
        out_dir=lambda wildcards: f"outputs/truvari/{wildcards.sample}.{wildcards.aligner}.{wildcards.caller}.{wildcards.svtype}",
    shell:
        "(truvari bench "
        "--base {input.truth_set} "
        "--comp {input.comp_file} "
        "--output {params.out_dir} "
        "--includebed {input.truth_set_regions} "
        "--pctseq 0 "
        "--sizemin 50 "
        "--sizemax 1_000_000 "
        "--passonly "
        "&& cat {params.out_dir}/summary.json > {output} "
        ")2> {log}"
