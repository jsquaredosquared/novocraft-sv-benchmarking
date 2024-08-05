configfile: "../../config/config.yaml"


SAMPLES = config["samples"]
ALIGNERS = config["aligners"]
CALLERS = config["callers"]
SVTYPES = config["truth_set_svtypes"]


rule all:
    input:
        "test.txt"


rule generate_plots:
    input:
        collect("../../outputs/truvari/{sample}.{aligner}.{caller}.{svtype}.truvari-bench.json", sample=SAMPLES, aligner=ALIGNERS, caller=CALLERS, svtype=SVTYPES)
    output:
        "test.txt"
    conda:
        "../envs/data_analysis.yaml"
    notebook:
        "../notebooks/generate_plots.py.ipynb"