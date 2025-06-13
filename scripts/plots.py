"""Create Manhattan and QQ plots from EWAS results."""

import argparse
import os

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scipy.stats as stats

# Parse arguments
parser = argparse.ArgumentParser(description="Plot Manhattan and QQ plots from EWAS results")
parser.add_argument("--input-file", "-i", required=True, help="Annotated EWAS results file")
parser.add_argument("--out-dir", required=True, help="Directory to save plots")
parser.add_argument("--assoc", required=True, help="Association variable")
args = parser.parse_args()

# Load annotated EWAS results
df = pd.read_csv(args.input_file, sep=None, engine="python")

# Calculate -log10(p-value)
df["-log10(pvalue)"] = -np.log10(df["pvalue"])

# Manhattan plot
fig_manhattan = px.scatter(
    df,
    x="MAPINFO" if "MAPINFO" in df.columns else "position",
    y="-log10(pvalue)",
    color="CHR" if "CHR" in df.columns else None,
    hover_data=["CpG", "gene"] if "gene" in df.columns else ["CpG"],
    title=f"Manhattan Plot - {args.assoc}",
)

# QQ plot
df_sorted = df[["pvalue"]].dropna().sort_values("pvalue")
df_sorted["expected"] = -np.log10(
    stats.uniform.ppf((np.arange(1, len(df_sorted) + 1)) / (len(df_sorted) + 1))
)
df_sorted["observed"] = -np.log10(df_sorted["pvalue"].values)

lambda_gc = np.median(stats.chi2.isf(df_sorted["pvalue"], 1)) / 0.456

fig_qq = go.Figure()
fig_qq.add_trace(
    go.Scatter(x=df_sorted["expected"], y=df_sorted["observed"], mode="markers", name="Observed")
)
fig_qq.add_trace(
    go.Scatter(
        x=df_sorted["expected"],
        y=df_sorted["expected"],
        mode="lines",
        line=dict(dash="dash"),
        name="Expected",
    )
)
fig_qq.update_layout(
    title=f"QQ Plot - {args.assoc} (lambda={lambda_gc:.2f})",
    xaxis_title="Expected -log10(p)",
    yaxis_title="Observed -log10(p)",
)

# Save plots
os.makedirs(args.out_dir, exist_ok=True)
fig_manhattan.write_html(os.path.join(args.out_dir, f"manhattan_{args.assoc}.html"))
fig_qq.write_html(os.path.join(args.out_dir, f"qq_{args.assoc}.html"))

print(f"Saved Manhattan and QQ plots to {args.out_dir}")
