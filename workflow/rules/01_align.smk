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
        expand("../../resources/bam-files/{sample}.{aligner}.bam", sample=SAMPLES, aligner=ALIGNERS),


rule align_with_bwa_mem:
    input:
        get_input_fastqs
    output:
        "../../resources/bam-files/{sample}.bwa-mem.bam",
    log:
        "../../logs/align_{sample}_with_bwa-mem.log"
    threads: 28
    conda:
        "../envs/bwa_env.yaml"
    shell:
        "bwa-mem2 mem "
        "-t 24 "
        f"{REFERENCE} "
        "{input} | "
        "samtools sort -o {output} -@ 4 "
        "2> {log}"


rule align_with_novoalign:
    input:
        get_input_fastqs,
    output:
        "../../resources/bam-files/{sample}.novoalign.bam",
    log:
        "../../logs/align_{sample}_with_novoalign.log"
    params:
        ref_nix = f"{REFERENCE}.nix"
    conda:
        "../envs/bwa_env.yaml"
    threads: 20
    shell:
        f"{ALIGNERS["novoalign"]} "
        "-i 400,100 "
        f"-d {NOVOINDEX} "
        "-f {input} "
        "-o BAM | "
        "samtools sort -o {output} -@ 4 "
        "2> {log}"

