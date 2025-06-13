# make_bed.py - Python version of make_bed.R

import argparse
import pandas as pd
import os

# Parse arguments
parser = argparse.ArgumentParser(description="Create a BED-style file from EWAS results")
parser.add_argument("--results", required=True, help="Path to annotated EWAS results")
parser.add_argument("--out-dir", default="results/bed", help="Output directory")
parser.add_argument("--assoc", required=True, help="Name of the association variable")
parser.add_argument("--p-threshold", type=float, default=1e-5, help="P-value threshold for filtering CpGs")
args = parser.parse_args()

# Load data
df = pd.read_csv(args.results, sep=None, engine='python')

# Filter by p-value
df_filtered = df[df['pvalue'] < args.p_threshold].copy()

# Ensure required columns exist
required_cols = ['CHR', 'MAPINFO', 'CpG', 'pvalue']
missing = [col for col in required_cols if col not in df_filtered.columns]
if missing:
    raise ValueError(f"Missing required columns in EWAS results: {missing}")

# BED format: chrom, start, end, name, score (we use pvalue as score)
df_bed = pd.DataFrame({
    'chrom': df_filtered['CHR'].astype(str),
    'start': df_filtered['MAPINFO'] - 1,
    'end': df_filtered['MAPINFO'],
    'name': df_filtered['CpG'],
    'score': df_filtered['pvalue']
})

# Save to BED-like file
os.makedirs(args.out_dir, exist_ok=True)
out_path = os.path.join(args.out_dir, f"{args.assoc}_ewas_results.bed")
df_bed.to_csv(out_path, sep='\t', index=False, header=False)
print(f"Saved BED file to {out_path}")
