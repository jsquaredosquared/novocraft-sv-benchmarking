configfile: "../../config/config.yaml"


INPUT_PREFIX = config["input_prefix"]
OUTPUT_PREFIX = config["output_prefix"]
SAMPLES = config["samples"]
ALIGNERS = config["aligners"]
REFERENCE = config["reference"]
NOVOINDEX = config["novoindex"]


def get_input_fastqs(wildcards):
    return [f"{INPUT_PREFIX}/{fastq}" for fastq in SAMPLES[wildcards.sample]]


rule all:
    input:
        expand("../../resources/bam-files/{sample}.{aligner}_sorted.bam", sample=SAMPLES, aligner=ALIGNERS),


rule align_with_bwa_mem:
    input:
        get_input_fastqs
    output:
        temp("../../resources/bam-files/{sample}.bwa-mem.sam"),
    log:
        "../../logs/align_{sample}_with_bwa-mem.log"
    threads: 24
    container:
        "docker://dceoy/bwa-mem2:latest"
    shell:
        "bwa-mem2 mem "
        "-t 24 "
        f"{REFERENCE} "
        "{input} "
        "> {output} "
        "2> {log}"


rule align_with_novoalign:
    input:
        get_input_fastqs,
    output:
        temp("../../resources/bam-files/{sample}.novoalign.sam"),
    log:
        "../../logs/align_{sample}_with_novoalign.log"
    params:
        ref_nix = f"{REFERENCE}.nix"
    threads: 24
    shell:
        f"{ALIGNERS["novoalign"]} "
        "-i 400,100 "
        f"-d {NOVOINDEX} "
        "-f {input} "
        "-o SAM "
        "{output} "
        "2> {log}"

rule sort_sam_to_bam:
    input:
        "../../resources/bam-files/{sample}.{aligner}.sam"
    output:
        "../../resources/bam-files/{sample}.{aligner}_sorted.bam"
    container:
        "docker://staphb/samtools:latest"
    log:
        "../../logs/sort_{sample}_{aligner}_sam_to_bam.log"
    threads: 16
    shell:
        "samtools sort {input} -o {output} -@ 16 2> {log}"

