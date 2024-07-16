configfile: "../../config/config.yaml"

SAMPLES = config["samples"]
CALLERS = config["callers"]
REFERENCE = config["reference"]

rule all:
    input: 
        expand("../../outputs/{sample}.{aligner}.{caller}.vcf", sample=SAMPLES, aligner=ALIGNERS, caller=CALLERS)

rule call_svs_with_manta:
    input:
        "../../resources/bam-files/{sample}.{aligner}.bam"
    output:
        "../../outputs/{sample}.{aligner}.{caller}.vcf"
    conda:
        "manta_env.yaml"
    log:
        "../../logs/{sample}.{aligner}.{caller}.log"
    params:
        out_dir = "../../outputs/{caller}"
    shell:
        "configManta.py "
        "--bam {input} "
        f"--referenceFasta {REFERENCE} "
        "--runDir {params.out_dir} "
        "2> {log}"
    

