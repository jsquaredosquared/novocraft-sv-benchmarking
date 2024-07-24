

rule compare_to_truth_set:
    input:
        "../../outputs/{caller}/{sample}.{aligner}.{caller}.vcf"
    output:
        ""
    conda:
    log:
    params:
    shell:
        "truvari bench "
        "-b ??? "
        "--includebed ???"
        "-c {input} "
        "-o ??? "
        "--sizemin 50 "
        "--sizemax 1_000_000 "