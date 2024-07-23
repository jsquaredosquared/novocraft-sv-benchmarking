configfile: "../../config/config.yaml"


SAMPLES = config["samples"]
CALLERS = config["callers"]
ALIGNERS = config["aligners"]
REFERENCE = config["reference"]


rule all:
    input:
        expand(
            "../../outputs/dysgu/{sample}.{aligner}.dysgu.vcf",
            sample=SAMPLES,
            aligner=ALIGNERS
        )


rule configure_manta:
    input:
        "../../resources/alignment-files/{sample}.{aligner}_sorted.cram",
    output:
        "../../outputs/manta/{sample}/{aligner}/runWorkflow.py"
    conda:
        "../envs/manta_env.yaml"
    log:
        "../../logs/{sample}_{aligner}.configure_manta.log",
    params:
        out_dir = lambda wildcards: f"../../outputs/{wildcards.caller}/{wildcards.sample}/{wildcards.aligner}",
        regions_bed = "../../resources/manta/manta_main-contigs.bed.gz"
    shell:
        "configManta.py "
        "--bam {input} "
        f"--referenceFasta {REFERENCE} "
        "--callRegions {params.regions_bed} "
        "--runDir {params.out_dir} "
        "2> {log}"


rule run_manta:
    input:
        "../../outputs/manta/{sample}/{aligner}/runWorkflow.py"
    output:
        "../../outputs/manta/{sample}/{aligner}/diploidSV.vcf.gz",
        "../../outputs/manta/{sample}/{aligner}/candidateSV.vcf.gz",
        "../../outputs/manta/{sample}/{aligner}/candidateSmallIndels.vcf.gz",
    conda:
        "../envs/manta_env.yaml"
    log:
        "../../logs/{sample}_{aligner}.execute_manta.log"
    threads: 30
    shell:
        "{input} "
        "-j 30 "
        "-e jeffrey@novocraft.com"


rule run_dysgu:
    input:
        "../../resources/alignment-files/{sample}.{aligner}_sorted.cram"
    output:
        "../../outputs/dysgu/{sample}.{aligner}.dysgu.vcf"
    conda:
        "../envs/dysgu_env.yaml"
    log:
        "../../logs/{sample}.{aligner}.dysgu.log"
    threads: 30
    params:
        temp_dir = "../../outputs/dysgu/{sample}.{aligner}.temp"
    shell:
        "dysgu run --clean -p30 "
        f"{REFERENCE} "
        "{params.temp_dir} "
        "{input} "
        "> {output} "
        "2> {log}"