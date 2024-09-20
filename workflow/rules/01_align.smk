def get_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]


def get_reference(wildcards):
    return config["reference"]


rule generate_cram_files:
    input:
        expand(
            "resources/alignment-files/{sample}.{aligner}.{ext}",
            sample=config["samples"],
            aligner=config["aligners"],
            ext=[".cram", ".cram.crai"],
        ),


rule align_with_bwa_mem2:
    input:
        fastq = get_input_fastqs,
        reference = get_reference
    output:
        "resources/alignment-files/{sample}.bwa-mem2.cram",
    log:
        "logs/align_{sample}_with_bwa-mem2.log",
    conda:
        "../envs/alignment.yaml"
    threads: 32
    shell:
        "(bwa-mem2 mem "
        "-t 24 "
        "-R '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA' "
        "{input.reference} {input.fastq} "
        "| samtools sort -@ 4 -O bam -l 0 -T /tmp - "
        "| samtools view -@ 4 -T {input.reference} -C -o {output} - "
        ")2> {log}"


rule align_with_novoalign:
    input:
        fastq = get_input_fastqs,
        reference = get_reference,
    output:
        "resources/alignment-files/{sample}.novoalign.cram",
    log:
        "logs/align_{sample}_with_novoalign.log",
    conda:
        "../envs/alignment.yaml"
    threads: 48
    shell:
        "({config[aligners][novoalign]} -d {input.reference}.nix -c 40 "
        "-f {input.fastq} -o SAM "
        "'@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA' "
        "| samtools sort -@ 4 -O bam -l 0 -T /tmp - "
        "| samtools view -@ 4 -T {input.reference} -C -o {output} - "
        ")2> {log}"


rule index_cram_file:
    input:
        "resources/alignment-files/{sample}.{aligner}.cram",
    output:
        "resources/alignment-files/{sample}.{aligner}.cram.crai",
    log:
        "logs/index_{sample}.{aligner}_cram.log",
    conda:
        "../envs/alignment.yaml"
    threads: 8
    shell:
        "samtools index {input} -@ {threads}"


# http://www.htslib.org/workflow/cram.html
# TODO: Should you mark duplicates?
