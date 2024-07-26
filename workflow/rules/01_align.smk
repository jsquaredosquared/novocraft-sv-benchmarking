configfile: "../../config/config.yaml"


SAMPLES = config["samples"]
ALIGNERS = config["aligners"]
REFERENCE = config["reference"]


def get_input_fastqs(wildcards):
    return SAMPLES[wildcards.sample]


rule all:
    input:
        expand(
            "../../resources/alignment-files/{sample}.{aligner}.{ext}",
            sample=SAMPLES,
            aligner=ALIGNERS,
            ext=["cram", "cram.crai"]
        ),


rule align_with_bwa_mem2:
    input:
        get_input_fastqs
    output:
        "../../resources/alignment-files/{sample}.bwa-mem2.cram",
    log:
        "../../logs/align_{sample}_with_bwa-mem2.log",
    conda:
        "../envs/alignment_env.yaml"
    threads: 32
    shell:
        f"(bwa-mem2 mem "
        "-t 24 "
        "-R '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA'"
        f"{REFERENCE} {{input}} "
        "| samtools sort -@ 4 -O bam -l 0 -T /tmp - "
        f"| samtools view -@ 4 -T {REFERENCE} -C -o {{output}} - "
        ")2> {log}"


rule align_with_novoalign:
    input:
        get_input_fastqs
    output:
        "../../resources/alignment-files/{sample}.novoalign.cram",
    log:
        "../../logs/align_{sample}_with_novoalign.log",
    conda:
        "../envs/alignment_env.yaml"
    threads: 32
    shell:
        f"({ALIGNERS["novoalign"]} -d {REFERENCE}.nix "
        "-f {input} -o SAM "
        "'@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA' "
        "| samtools sort -@ 4 -O bam -l 0 -T /tmp - "
        f"| samtools view -@ 4 -T {REFERENCE} -C -o {{output}} - "
        ")2> {log}"
        

rule index_cram_file:
    input:
        "../../resources/alignment-files/{sample}.{aligner}.cram"
    output:
        "../../resources/alignment-files/{sample}.{aligner}.cram.crai"
    log:
        "../../logs/index_{sample}.{aligner}_cram.log"
    conda:
        "../envs/alignment_env.yaml"
    threads: 16
    shell:
        "samtools index {input} -@ {threads}"
