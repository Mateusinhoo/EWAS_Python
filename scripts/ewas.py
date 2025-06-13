# ewas.py - Python version of ewas.R

import argparse
import pandas as pd
import numpy as np
from statsmodels.api import OLS, add_constant
from concurrent.futures import ProcessPoolExecutor
import os

# Parse arguments
parser = argparse.ArgumentParser(description="Run EWAS analysis")
parser.add_argument("--pheno", required=True, help="Phenotype file")
parser.add_argument("--methyl", required=True, help="Methylation matrix")
parser.add_argument("--assoc", required=True, help="Variable to associate")
parser.add_argument("--chunk-size", type=int, default=1000)
parser.add_argument("--out-dir", default="results/")
parser.add_argument("--out-type", default=".csv")
parser.add_argument("--workers", type=int, default=4)
args = parser.parse_args()

# Load data
pheno = pd.read_csv(args.pheno)
mvals = pd.read_csv(args.methyl, index_col=0)

assert args.assoc in pheno.columns, f"Association variable {args.assoc} not found"

# Merge data
common_samples = pheno["sample_id"].isin(mvals.columns)
pheno = pheno[common_samples].set_index("sample_id")
mvals = mvals[pheno.index]

# Define analysis function
def analyze_cpg(cpg_id):
    y = mvals.loc[cpg_id]
    X = add_constant(pheno[[args.assoc]].astype(float))
    try:
        model = OLS(y, X).fit()
        return {
            "CpG": cpg_id,
            "beta": model.params[args.assoc],
            "pvalue": model.pvalues[args.assoc]
        }
    except Exception:
        return {"CpG": cpg_id, "beta": np.nan, "pvalue": np.nan}

# Run analysis in chunks
results = []
cpg_list = mvals.index.tolist()
with ProcessPoolExecutor(max_workers=args.workers) as executor:
    for res in executor.map(analyze_cpg, cpg_list):
        results.append(res)

# Save output
out_df = pd.DataFrame(results)
os.makedirs(args.out_dir, exist_ok=True)
out_path = os.path.join(args.out_dir, f"ewas_results{args.out_type}")
out_df.to_csv(out_path, index=False)
print(f"Results saved to {out_path}")
