configfile: "../../config/config.yaml"


INPUT_PREFIX = config["prefix"]
OUTPUT_PREFIX = config["results"]
SAMPLES = config["fastq"]
ALIGNERS = config["aligners"]
REFERENCE = config["reference"]


def get_input_fastqs(wildcards):
    return [f"{INPUT_PREFIX}/{fastq}" for fastq in SAMPLES[wildcards.sample]]


rule all:
    input:
        expand("{sample}.{aligner}.bam", sample=SAMPLES, aligner=ALIGNERS),


rule align_with_bwa_mem:
    input:
        get_input_fastqs,
    output:
        "{sample}.bwa-mem.bam",
    log:
        "../../logs/align_{sample}_with_bwa-mem.log"
    params:
        reference=REFERENCE,
        bwa_mem=ALIGNERS["bwa-mem"],
        samtools="docker run staphb/samtools"
    threads:
    shell:
        "({bwa_mem} mem "
        "-t ?? "
        "{reference} "
        "{input} | "
        "{samtools} sort -o {output} -@??) "
        "2> {log}"


rule align_with_novoalign:
    input:
        get_input_fastqs,
    output:
        "{sample}.novoalign.bam",
    log:
        "../../logs/align_{sample}_with_novoalign.log"
    params:
        novoalign=ALIGNERS["novoalign"],
    shell:
        "{novoalign} "
        "-i 400,100 "
        "-d ???index??? "
        "-f {input} "
        "-o BAM "
        "> {output} "
        "2> {log}"

