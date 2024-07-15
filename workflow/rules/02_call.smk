configfile: "../../config/config.yaml"

OUTPUT_PREFIX = "../../outputs/manta"
SAMPLES = config["samples"]
CALLERS = config["callers"]

rule all:
    input: 
        expand("../../outputs/{sample}.{caller}.vcf", sample=SAMPLES, caller=CALLERS)

rule call_svs_with_manta:
    input:
    output:
    conda:
    shell:
    

