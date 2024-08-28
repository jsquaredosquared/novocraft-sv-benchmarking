rule generate_plots:
    input:
        collect(
            "outputs/truvari/{sample}.{aligner}.{caller}.{svtype}.truvari-bench.json",
            sample=config["samples"],
            aligner=config["aligners"],
            caller=config["callers"],
            svtype=config["truth_set_svtypes"],
        ),
    output:
        "results/final_plot.svg",
    conda:
        "workflow/envs/data_analysis.yaml"
    script:
        "scripts/generate_plots.py"
