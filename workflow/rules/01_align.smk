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
            "../../resources/bam-files/{sample}.{aligner}_sorted.bam",
            sample=SAMPLES,
            aligner=ALIGNERS,
        ),


rule align:
    input:
        get_input_fastqs,
    output:
        temp("../../resources/bam-files/{sample}.{aligner}.sam"),
    log:
        "../../logs/align_{sample}_with_{aligner}.log",
    threads: 24
    run:
        match wildcards.aligner:
            case "bwa-mem":
                shell(
                    f"bwa-mem2 mem -t 24 {REFERENCE} {{input}} > {{output}} 2> {{log}}"
                )
            case "novoalign":
                shell(
                    f"{ALIGNERS["novoalign"]} -i 400,100 -d {REFERENCE}.nix -f {{input}} -o SAM '@RG\tID:V3\tSM:NA12878\tPL:ILLUMINA\tLB:sv' > {{output}} 2> {{log}}"
                )


rule sort_sam_to_cram:
    input:
        "../../resources/bam-files/{sample}.{aligner}.sam",
    output:
        "../../resources/bam-files/{sample}.{aligner}_sorted.cram",
    log:
        "../../logs/sort_{sample}_{aligner}_sam_to_cram.log",
    threads: 16
    shell:
        "samtools sort {input} -o {output} -@ 16 2> {log}"
