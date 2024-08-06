configfile: "../../config/config.yaml"


rule all:
    input:
        "test.txt",


rule generate_plots:
    input:
        collect(
            "../../outputs/truvari/{sample}.{aligner}.{caller}.{svtype}.truvari-bench.json",
            sample=SAMPLES,
            aligner=ALIGNERS,
            caller=CALLERS,
            svtype=SVTYPES,
        ),
    output:
        "../../results/final_plot.svg",
    conda:
        "../envs/data_analysis.yaml"
    notebook:
        "../notebooks/generate_plots.py.ipynb"
