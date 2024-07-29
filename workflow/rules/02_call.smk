configfile: "../../config/config.yaml"


SAMPLES = config["samples"]
CALLERS = config["callers"]
# ALIGNERS = config["aligners"]
ALIGNERS = ["bwa-mem2"]  # Couldn't wait for novoalign to finish
REFERENCE = config["reference"]


rule all:
    input:
        expand(
            "../../outputs/{caller}/{sample}.{aligner}.{caller}.vcf",
            sample=SAMPLES,
            aligner=ALIGNERS,
            caller=CALLERS,
        ),


rule configure_manta:
    input:
        "../../resources/alignment-files/{sample}.{aligner}_sorted.cram",
    output:
        "../../outputs/manta/{sample}/{aligner}/runWorkflow.py",
    conda:
        "../envs/manta_env.yaml"
    log:
        "../../logs/{sample}_{aligner}.configure_manta.log",
    params:
        out_dir=lambda wildcards: f"../../outputs/manta/{wildcards.sample}/{wildcards.aligner}",
        regions_bed="../../resources/manta/manta_main-contigs.bed.gz",
    shell:
        "configManta.py "
        "--bam {input} "
        f"--referenceFasta {REFERENCE} "
        "--callRegions {params.regions_bed} "
        "--runDir {params.out_dir} "
        "2> {log}"


# TODO Use script to reformat inversions into single fields
rule run_manta:
    input:
        "../../outputs/manta/{sample}/{aligner}/runWorkflow.py",
    output:
        "../../outputs/manta/{sample}.{aligner}.manta.vcf"
    conda:
        "../envs/manta_env.yaml"
    log:
        "../../logs/{sample}_{aligner}.execute_manta.log",
    threads: 16
    shell:
        "(python {input} "
        "-j {threads} "
        "-e jeffrey@novocraft.com "
        "&& gunzip --stdout ../../outputs/manta/{wildcards.sample}/{wildcards.aligner}/diploidSV.vcf.gz > {output} "
        ")2> {log}"


rule run_dysgu:
    input:
        "../../resources/alignment-files/{sample}.{aligner}.cram",
    output:
        "../../outputs/dysgu/{sample}.{aligner}.dysgu.vcf",
    conda:
        "../envs/dysgu_env.yaml"
    log:
        "../../logs/{sample}.{aligner}.dysgu.log",
    threads: 16
    params:
        temp_dir="../../outputs/dysgu/{sample}.{aligner}.temp",
    shell:
        "dysgu run --clean -p {threads} "
        f"{REFERENCE} "
        "{params.temp_dir} "
        "{input} "
        "> {output} "
        "2> {log}"


# No multithreading (https://github.com/dellytools/delly/issues/268#issuecomment-975385454).
rule run_delly:
    input:
        "../../resources/alignment-files/{sample}.{aligner}.cram",
    output:
        "../../outputs/delly/{sample}.{aligner}.delly.vcf",
    conda:
        "../envs/delly_env.yaml"
    log:
        "../../logs/{sample}.{aligner}.delly.log",
    threads: 1
    params:
        exclude_regions="../../resources/delly/human.hg19.excl.tsv",
    shell:
        f"delly call -g {REFERENCE} "
        "-x {params.exclude_regions} "
        "{input} > {output} "
        "2> {log}"


rule run_smoove:
    input:
        "../../resources/alignment-files/{sample}.{aligner}.cram",
    output:
        "../../outputs/lumpy/{sample}.{aligner}.lumpy.vcf",
    conda:
        "../envs/smoove_env.yaml"
    log:
        "../../logs/{sample}.{aligner}.lumpy.log",
    params:
        exclude_regions="../../resources/lumpy/ceph18.b37.lumpy.exclude.2014-01-15.bed",
    threads: 16
    shell:
        "(smoove call "
        "--name {wildcards.sample}.{wildcards.aligner} "
        "--genotype "
        "--removepr "
        "--duphold "
        "--exclude {params.exclude_regions} "
        f"--fasta {REFERENCE} "
        "-p {threads} "
        "--outdir ../../outputs/lumpy "
        "{input} "
        "&& gunzip --stdout ../../outputs/lumpy/{wildcards.sample}-smoove.genotyped.vcf.gz > {output} "
        ")2> {log}"


# TODO: Should you use the filter script provided by the devs?
rule run_wham:
    input:
        "../../resources/alignment-files/{sample}.{aligner}.cram",
    output:
        "../../outputs/wham/{sample}.{aligner}.wham.vcf",
    conda:
        "../envs/wham_env.yaml"
    log:
        "../../logs/{sample}.{aligner}.wham.log",
    params:
        exclude="GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605",
    threads: 10
    shell:
        "whamg -f {input} "
        f"-a {REFERENCE} "
        "-e {params.exclude} "
        "-x {threads} "
        "> {output} "
        "2> {log} "


rule run_tiddit:
    input:
        "../../resources/alignment-files/{sample}.{aligner}.cram",
    output:
        "../../outputs/tiddit/{sample}.{aligner}.tiddit.vcf",
    conda:
        "../envs/tiddit_env.yaml"
    log:
        "../../logs/{sample}.{aligner}.tiddit.log"
    threads: 1
    shell:
        "tiddit --sv --bam {input} "
        "-o ../../outputs/tiddit/{wildcards.sample}.{wildcards.aligner}.tiddit "
        f"--ref {REFERENCE} "
        "2> {log}"