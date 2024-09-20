import json
import re
from statistics import harmonic_mean

import altair as alt
import numpy as np
import polars as pl
import polars.selectors as cs
from snakemake.script import snakemake

alt.data_transformers.enable("vegafusion")

# region PREPARE DATA
records = []

for bench_file in snakemake.input:
    record = {}

    data = json.load(open(bench_file))

    record["f1"] = data["f1"]
    record["recall"] = data["recall"]
    record["precision"] = data["precision"]

    record["sample"], record["aligner"], record["caller"], record["svtype"] = re.search(
        r"(\w+).([\w\-_]+).(\w+).(ALL|DEL|INS).truvari-bench", bench_file
    ).groups()

    records.append(record)

benchmark_df = pl.DataFrame(records)

benchmark_all = benchmark_df.filter(pl.col("svtype") == "ALL")


final_table = benchmark_all.select(
    ["sample", "aligner", "caller", "svtype", "f1", "recall", "precision"]
).sort("f1", descending=True)
# endregion


# region PLOT DATA

# region Plot base precision-recall plot

lines_of_change = (
    alt.Chart(benchmark_all)
    .mark_line()
    .encode(
        x="recall",
        y="precision",
        color=alt.Color("caller", legend=None).scale(scheme="set1"),
    )
)

precision_recall_points = (
    alt.Chart(benchmark_all)
    .mark_point()
    .encode(
        x=alt.X("recall:Q").scale(domain=(0.1, 0.55)),
        y=alt.Y("precision:Q").scale(domain=(0.5, 1.0)),
        shape="aligner",
        fill=alt.Fill("caller").scale(scheme="set1"),
        stroke=alt.Stroke("caller").scale(scheme="set1"),
        tooltip=["aligner", "caller", "recall", "precision", "f1"],
    )
    .properties(width=500, height=500)
)

precision_recall_plot = lines_of_change + precision_recall_points

# endregion

# region Create F1 contours

zero_to_one = np.arange(0, 1, 0.001)

# Generate coordinates for each F1 contour.
# Recall and F1 are known, so we need to determine precision by first
# rearranging the F1 formula.

f1_scores = np.arange(0.1, 1.01, 0.1)

f1_contours = (
    pl.DataFrame(
        {
            "f1": np.array([[f1] * 1000 for f1 in f1_scores]).flatten(),
            "recall": list(zero_to_one) * 10,
        }
    )
    .filter(2 * pl.col("recall") - pl.col("f1") != 0)
    .with_columns(
        (pl.col("f1") * pl.col("recall") / (2 * pl.col("recall") - pl.col("f1"))).alias(
            "precision"
        )
    )
    .filter(pl.col("precision") > 0)
)

f1_contour_plot = alt.layer(
    *[
        (
            alt.Chart(f1_contours.filter(pl.col("f1") == f1))
            .mark_line()
            .encode(
                x=alt.X("recall"),
                y=alt.Y("precision").scale(domain=(0.5, 1)),
                color=alt.value("#808080"),
                tooltip=["f1"],
            )
        )
        for f1 in f1_scores
    ]
)

# endregion

# region Create contour labels

f1_labels = np.arange(0.1, 1.01, 0.1).round(1)

# For final zoomed-in image, display labels at precision (y) = 0.85
precision_labels = [0.85] * len(f1_labels)

# Need to find recall (x) coordinate for label.
label_coords = (
    pl.DataFrame({"f1": f1_labels, "precision": precision_labels})
    .with_columns(
        (
            pl.col("f1")
            * pl.col("precision")
            / (2 * pl.col("precision") - pl.col("f1"))
        ).alias("recall")
    )
    .with_columns(
        (
            (-pl.col("f1") ** 2 / (2 * pl.col("recall") - pl.col("f1")) ** 2)
            .arctan()
            .degrees()
            .abs()
            .round()
        ).alias("angle")
    )
    .with_columns(
        (
            pl.col("f1").map_elements(lambda f1: f"F1={f1}", return_dtype=pl.String)
        ).alias("f1_text")
    )
)

label_chart = (
    alt.Chart(label_coords)
    .mark_text(xOffset=8)
    .encode(x="recall", y="precision", text="f1_text", angle=alt.value(80))
)

# endregion

# region Compose final plot

final_plot = (
    alt.layer(f1_contour_plot, label_chart, precision_recall_plot)
    .properties(title="HG002 SV benchmarking results")
    .configure_line(strokeWidth=0.8, color="#000000", opacity=0.7)
    .configure_point(size=75)
    .configure_title()
).interactive()


final_plot.save(snakemake.output)

# endregion
