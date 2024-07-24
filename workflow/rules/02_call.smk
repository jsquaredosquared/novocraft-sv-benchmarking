configfile: "../../config/config.yaml"


SAMPLES = config["samples"]
CALLERS = config["callers"]
ALIGNERS = config["aligners"]
REFERENCE = config["reference"]


rule all:
    input:
        expand(
            "../../outputs/{caller}/{sample}.{aligner}.{caller}.vcf",
            sample=SAMPLES,
            aligner=ALIGNERS,
            caller=CALLERS
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
    shell:output:
        "test.txt"
        "configManta.py "
        "--bam {input} "
        f"--referenceFasta {REFERENCE} "
        "--callRegions {params.regions_bed} "
        "--runDir {params.out_dir} "
        "2> {log}"


# TODO Figure out how to run with Python2
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
    container:
        "docker://kcleal/dysgu"
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


# No multithreading (https://github.com/dellytools/delly/issues/268#issuecomment-975385454).
rule run_delly:
    input:
        "../../resources/alignment-files/{sample}.{aligner}_sorted.cram"
    output:
        "../../outputs/delly/{sample}.{aligner}.delly.vcf"
    conda:
        "../envs/delly_env.yaml"
    log:
        "../../logs/{sample}.{aligner}.delly.log"
    threads:
    params:
        exclude_regions = "../../resources/delly/human.hg19.excl.tsv"
    shell:
        f"delly call -g {REFERENCE} "
        "-x {params.exclude_regions} "
        "{input} > {output} "
        "2> {log}"


rule run_smoove:
    input:
        "../../resources/alignment-files/{sample}.{aligner}_sorted.cram"
    output:
        "../../outputs/lumpy/{sample}.{aligner}.lumpy.vcf"
    conda:
        "../envs/smoove_env.yaml"
    params:
        exclude_regions = "../../resources/lumpy/ceph18.b37.lumpy.exclude.2014-01-15.bed"
    threads: 30
    shell:
        "smoove call "
        "-x ??? "
        "--genotype "
        "--name ??? "
        "--exclude {exclude_regions} "
        f"--fasta {REFERENCE} "
        "-p {threads} "
        "--outdir ??? "
        "{input} "