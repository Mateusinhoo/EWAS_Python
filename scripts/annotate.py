"""Annotate EWAS results with genomic information.

This script mimics the functionality of the original ``annotation.R`` script. It
reads an EWAS results file, merges the results with CpG annotations and writes a
new annotated result file.
"""

import argparse
import os

import pandas as pd

# Parse arguments
parser = argparse.ArgumentParser(description="Add annotation data to EWAS results")
parser.add_argument("--input-file", "-i", required=True, help="EWAS result file (CSV or TSV)")
parser.add_argument("--out-dir", required=True, help="Directory to save annotated results")
parser.add_argument(
    "--stratified", choices=["yes", "no"], default="no", help="Was the analysis stratified"
)
parser.add_argument("--assoc", required=True, help="Association variable name")
parser.add_argument("--out-type", default=".csv", help="Output file extension")
args = parser.parse_args()

# Load EWAS results
df = pd.read_csv(args.input_file, sep=None, engine="python")

# Load annotation files (assumed to be in annotation_files directory)
anno_cpg = pd.read_csv("annotation_files/EPIC_hg38.tsv.gz", sep="\t")
anno_snp = pd.read_csv("annotation_files/EPIC_snp_key.tsv.gz", sep="\t")

# Merge EWAS results with CpG annotations
annotated = df.merge(anno_cpg, on="CpG", how="left")

# Save annotated output
os.makedirs(args.out_dir, exist_ok=True)
out_path = os.path.join(args.out_dir, f"{args.assoc}_ewas_annotated_results{args.out_type}")
annotated.to_csv(out_path, index=False)
print(f"Annotated results saved to {out_path}")
