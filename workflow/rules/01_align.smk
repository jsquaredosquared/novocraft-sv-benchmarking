configfile: "../../config/config.yaml"


INPUT_PREFIX = config["input_prefix"]
OUTPUT_PREFIX = config["output_prefix"]
SAMPLES = config["samples"]
ALIGNERS = config["aligners"]
REFERENCE = config["reference"]


def get_input_fastqs(wildcards):
    return [f"{INPUT_PREFIX}/{fastq}" for fastq in SAMPLES[wildcards.sample]]


rule all:
    input:
        expand(
            "../../resources/alignment-files/{sample}.{aligner}_sorted.bam",
            sample=SAMPLES,
            aligner=ALIGNERS,
        ),


rule align_with_bwa_mem:
    input:
        get_input_fastqs,
    output:
        temp("../../resources/alignment-files/{sample}.{aligner}.sam"),
    log:
        "../../logs/align_{sample}_with_{aligner}.log",
    conda:
        "../envs/alignment_env.yaml"
    threads: 24
    shell:
        f"bwa-mem2 mem -t 24 {REFERENCE} {{input}} > {{output}} 2> {{log}}"


rule align_with_novoalign:
    input:
        get_input_fastqs,
    output:
        temp("../../resources/alignment-files/{sample}.{aligner}.sam"),
    log:
        "../../logs/align_{sample}_with_{aligner}.log",
    threads: 24
    run:
        f"{ALIGNERS["novoalign"]} -i 400,100 -d {REFERENCE}.nix -f {{input}} -o SAM '@RG\tID:V3\tSM:NA12878\tPL:ILLUMINA\tLB:sv' > {{output}} 2> {{log}}"

rule sort_sam_to_cram:
    input:
        "../../resources/alignment-files/{sample}.{aligner}.sam",
    output:
        "../../resources/alignment-files/{sample}.{aligner}_sorted.cram",
    log:
        "../../logs/sort_{sample}_{aligner}_sam_to_cram.log",
    threads: 16
    conda:
        "../envs/alignment_env.yaml"
    shell:
        "samtools sort {input} -o {output} -@ 16 2> {log}"
        

rule index_cram_file:
    input:
        "../../resources/alignment-files/{sample}.{aligner}_sorted.cram"
    output:
        "../../resources/alignment-files/{sample}.{aligner}_sorted.cram.crai"
    log:
        "../../logs/index_{sample}.{aligner}_cram.log"
    threads: 16
    conda:
        "../envs/alignment_env.yaml"
    shell:
        "samtools index {output} -@ 16"
