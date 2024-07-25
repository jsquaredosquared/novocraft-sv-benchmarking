configfile: "../../config/config.yaml"


TRUTH_SETS = config["truth_sets"]
REFERENCE = config["reference"]


rule compare_to_truth_set:
    input:
        "../../outputs/{caller}/{sample}.{aligner}.{caller}.vcf"
    output:
        "../../outputs/truvari/{sample}.{aligner}.{caller}/summary.json"
    conda:
        "../envs/truvari_env.yaml"
        # Check Truvari installation guide for caveats.
    log:
        "../../logs/{sample}.{aligner}.{caller}.truvari.log"
    params:
        truth_set = lambda wildcards: TRUTH_SETS[wildcards.sample],
        out_dir = lambda wildcards: f"../../outputs/truvari/{wildcards.sample}.{wildcards.aligner}.{wildcards.caller}"
    shell:
        "truvari bench "
        "--base {params.truth_set} "
        "--comp {input} "
        "--output {params.out_dir} "
        "--includebed ../../../../resources/sv-benchmark/HG002/HG002_SVs_Tier1_v0.6.bed "
        "--pctseq 0 "
        "--sizemin 50 "
        "--sizemax 1_000_000 "
        "--passonly "
        "2> {log}"