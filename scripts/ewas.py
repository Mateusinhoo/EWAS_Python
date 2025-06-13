"""Run a basic epigenome-wide association study (EWAS).

This script mirrors the functionality of ``ewas.R`` using ``statsmodels`` to
perform ordinary least squares regression for each CpG.
"""
import argparse
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from statsmodels.api import OLS, add_constant


def main() -> None:
    """Entry point for command line execution."""

    parser = argparse.ArgumentParser(description="Run EWAS analysis")
    parser.add_argument("--pheno", required=True, help="Phenotype file")
    parser.add_argument("--methyl", required=True, help="Methylation matrix")
    parser.add_argument("--assoc", required=True, help="Variable to associate")
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--out-dir", default="results/")
    parser.add_argument("--out-type", default=".csv")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--sample-id-col",
        default="sampleID",
        help="Column name for sample identifiers in the phenotype file",
    )
    args = parser.parse_args()

    # Load data
    pheno = pd.read_csv(args.pheno)
    mvals = pd.read_csv(args.methyl, index_col=0)

    if args.assoc not in pheno.columns:
        raise ValueError(f"Association variable {args.assoc} not found")

    # Validate sample ID column
    if args.sample_id_col not in pheno.columns:
        raise ValueError(
            f"Phenotype file must contain a '{args.sample_id_col}' column"
        )

    sample_ids = pheno[args.sample_id_col].astype(str)
    # Determine orientation of methylation matrix
    cols_match = sample_ids.isin(mvals.columns.astype(str)).all()
    rows_match = sample_ids.isin(mvals.index.astype(str)).all()
    if not cols_match and rows_match:
        # Samples appear on rows -> transpose
        mvals = mvals.T
        print("Transposed methylation matrix to match sample orientation")
    elif not cols_match and not rows_match:
        raise ValueError("Sample IDs do not match methylation matrix dimensions")

    # Merge data
    pheno = pheno.set_index(args.sample_id_col)
    mvals = mvals.loc[:, pheno.index]

    design_matrix = add_constant(pheno[[args.assoc]].astype(float))

    # Define analysis function
    def analyze_cpg(cpg_id):
        """Return association statistics for a single CpG."""

        y = mvals.loc[cpg_id]
        try:
            model = OLS(y, design_matrix).fit()
            return {
                "CpG": cpg_id,
                "beta": model.params[args.assoc],
                "pvalue": model.pvalues[args.assoc],
            }
        except Exception:
            return {"CpG": cpg_id, "beta": np.nan, "pvalue": np.nan}

    # Run analysis in chunks
    results = []
    cpg_list = mvals.index.tolist()
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for i in range(0, len(cpg_list), args.chunk_size):
            batch = cpg_list[i : i + args.chunk_size]
            for res in executor.map(analyze_cpg, batch):
                results.append(res)

    # Save output
    out_df = pd.DataFrame(results)
    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, f"ewas_results{args.out_type}")
    out_df.to_csv(out_path, index=False)
    print(f"Results saved to {out_path}")


if __name__ == "__main__":
    main()

