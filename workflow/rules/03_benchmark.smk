configfile: "../../config/config.yaml"


SAMPLES = config["samples"]
TRUTH_SETS = config["truth_sets"]
REFERENCE = config["reference"]
ALIGNERS = config["aligners"]
CALLERS = config["callers"]
SVTYPES = config["truth_set_svtypes"]


rule all:
    input:
        expand(
            "../../outputs/truvari/{sample}.{aligner}.{caller}.{svtype}.truvari-bench.json",
            sample=SAMPLES,
            aligner=ALIGNERS,
            caller=CALLERS,
            svtype=SVTYPES,
        ),


def get_filter_expression(wildcards):
    match wildcards.svtype.split("-"):
        case [svtype]:
            return "True" if svtype == "ALL" else f'INFO["SVTYPE"] == "{svtype}"'
        case [min_len, max_len]:
            return f'{int(min_len)}<=INFO["SVLEN"]<={int(max_len)}'


for sample in SAMPLES:

    rule split_truth_set:
        input:
            truth_set=TRUTH_SETS[sample],
        output:
            [TRUTH_SETS[sample].replace("ALL", svtype) for svtype in SVTYPES],
        conda:
            "../envs/vembrane.yaml"
        shell:
            f"xonsh split-truth-set-by-svtype.xsh {{input}} {','.join(SVTYPES)}"


rule split_vcf_into_svtype:
    input:
        "../../outputs/{caller}/{sample}.{aligner}.{caller}.vcf.gz",
    output:
        vcf_file=temp(
            "../../outputs/{caller}/{sample}.{aligner}.{caller}.{svtype}.vcf.gz"
        ),
        index_file=temp(
            "../../outputs/{caller}/{sample}.{aligner}.{caller}.{svtype}.vcf.gz.tbi"
        ),
    conda:
        "../envs/vembrane.yaml"
    log:
        "../../logs/{sample}.{aligner}.{caller}.split_vcf_{svtype}.log",
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
        comp_file="../../outputs/{caller}/{sample}.{aligner}.{caller}.{svtype}.vcf.gz",
        index_file="../../outputs/{caller}/{sample}.{aligner}.{caller}.{svtype}.vcf.gz.tbi",
        truth_set=lambda wildcards: TRUTH_SETS[f"{wildcards.sample}"].replace("ALL", f"{wildcards.svtype}"),
    output:
        "../../outputs/truvari/{sample}.{aligner}.{caller}.{svtype}.truvari-bench.json",
    conda:
        "../envs/truvari.yaml"
    log:
        "../../logs/{sample}.{aligner}.{caller}.{svtype}.truvari-bench.log",
    params:
        out_dir=lambda wildcards: f"../../outputs/truvari/{wildcards.sample}.{wildcards.aligner}.{wildcards.caller}.{wildcards.svtype}",
        regions="../../resources/sv-benchmarks/HG002/HG002_SVs_Tier1_v0.6.bed",
    shell:
        "(truvari bench "
        "--base {input.truth_set} "
        "--comp {input.comp_file} "
        "--output {params.out_dir} "
        "--includebed {params.regions} "
        "--pctseq 0 "
        "--sizemin 50 "
        "--sizemax 1_000_000 "
        "--passonly "
        "&& cat {params.out_dir}/summary.json > {output} "
        ")2> {log}"
