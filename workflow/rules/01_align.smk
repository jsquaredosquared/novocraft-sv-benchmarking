configfile: "../../config/config.yaml"


SAMPLES = config["samples"]
ALIGNERS = config["aligners"]
REFERENCE = config["reference"]


def get_input_fastqs(wildcards):
    return SAMPLES[wildcards.sample]


rule all:
    input:
        expand(
            "../../resources/alignment-files/{sample}.{aligner}_sorted.cram",
            sample=SAMPLES,
            aligner=ALIGNERS,
        ),


rule align_with_bwa_mem2:
    input:
        get_input_fastqs
    output:
        temp("../../resources/alignment-files/{sample}.bwa-mem2.sam"),
    log:
        "../../logs/align_{sample}_with_bwa-mem2.log",
    conda:
        "../envs/alignment_env.yaml"
    threads: 30
    shell:
        f"bwa-mem2 mem -t 30 {REFERENCE} {{input}} > {{output}} 2> {{log}}"


rule align_with_novoalign:
    input:
        get_input_fastqs
    output:
        temp("../../resources/alignment-files/{sample}.novoalign.sam"),
    log:
        "../../logs/align_{sample}_with_novoalign.log",
    threads: 30
    shell:
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
        "samtools index {input} -@ 16"
