configfile: "../../config/config.yaml"


TRUTH_SETS = config["truth_sets"]


rule compare_to_truth_set:
    input:
        "../../outputs/{caller}/{sample}.{aligner}.{caller}.vcf"
    output:
        "../../outputs/truvari/{sample}.{aligner}.{caller}/summary.json"
    conda:
    log:
    params:
        truth_set = lambda wildcards: TRUTH_SETS[wildcards.sample],
        out_dir = lambda wildcards: f"../../outputs/truvari/{wildcards.sample}.{wildcards.aligner}.{wildcards.caller}"
    shell:
        "truvari bench "
        "-b {params.truth_set} "
        "--includebed ../../../../resources/hg002-giab-benchmark/HG002_SVs_Tier1_v0.6.bed "
        "-c {input} "
        "-o {params.out_dir} "
        "--sizemin 50 "
        "--sizemax 1_000_000 "