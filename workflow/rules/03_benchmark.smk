configfile: "../../config/config.yaml"

SAMPLES = config["samples"]
TRUTH_SETS = config["truth_sets"]
REFERENCE = config["reference"]
# ALIGNERS = config["aligners"]
ALIGNERS = ["bwa-mem2", "novoalign"]
# CALLERS = config["callers"]
CALLERS = ["dysgu"]


rule all:
    input:
        expand("../../outputs/truvari/{sample}.{aligner}.{caller}.truvari-bench.json", sample=SAMPLES, aligner=ALIGNERS, caller=CALLERS)


rule compare_to_truth_set:
    input:
        "../../outputs/{caller}/{sample}.{aligner}.{caller}.vcf.gz",
    output:
        "../../outputs/truvari/{sample}.{aligner}.{caller}.truvari-bench.json",
    conda:
        "../envs/truvari_env.yaml",
    log:
        "../../logs/{sample}.{aligner}.{caller}.truvari.log",
    params:
        truth_set=lambda wildcards: TRUTH_SETS[wildcards.sample],
        out_dir=lambda wildcards: f"../../outputs/truvari/{wildcards.sample}.{wildcards.aligner}.{wildcards.caller}/",
        regions="../../resources/sv-benchmarks/HG002/HG002_SVs_Tier1_v0.6.bed",
    shell:
        "(truvari bench "
        "--base {params.truth_set} "
        "--comp {input} "
        "--output {params.out_dir} "
        "--includebed {params.regions} "
        "--pctseq 0 "
        "--sizemin 50 "
        "--sizemax 1_000_000 "
        "--passonly "
        "&& cat {params.out_dir}/summary.json > {output} "
        ")2> {log}"
