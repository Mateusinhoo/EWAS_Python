"""Vectorized EWAS using NumPy operations."""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.api import add_constant


def run_vectorized(
    design: np.ndarray, y: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return betas and p-values for multiple CpGs simultaneously."""
    # design: n_samples x 2, with intercept and variable
    x_tx_inv = np.linalg.inv(design.T @ design)
    x_tx_inv_x_t = x_tx_inv @ design.T
    df_resid = design.shape[0] - design.shape[1]
    var_factor = x_tx_inv[1, 1]

    betas = x_tx_inv_x_t @ y
    residuals = y - design @ betas
    mse = np.sum(residuals ** 2, axis=0) / df_resid
    se_beta = np.sqrt(mse * var_factor)
    t_stats = betas[1] / se_beta
    p_values = 2 * stats.t.sf(np.abs(t_stats), df_resid)
    return betas[1], p_values


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run EWAS analysis using vectorized regression"
    )
    parser.add_argument("--pheno", required=True, help="Phenotype file")
    parser.add_argument("--methyl", required=True, help="Methylation matrix")
    parser.add_argument("--assoc", required=True, help="Variable to associate")
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--out-dir", default="results/")
    parser.add_argument("--out-type", default=".csv")
    parser.add_argument(
        "--sample-id-col",
        default="sampleID",
        help="Sample ID column in phenotype file",
    )
    args = parser.parse_args()

    pheno = pd.read_csv(args.pheno)
    mvals = pd.read_csv(args.methyl, index_col=0)

    if args.assoc not in pheno.columns:
        raise ValueError(f"Association variable {args.assoc} not found")
    if args.sample_id_col not in pheno.columns:
        raise ValueError(
            f"Phenotype file must contain a '{args.sample_id_col}' column"
        )

    sample_ids = pheno[args.sample_id_col].astype(str)
    cols_match = sample_ids.isin(mvals.columns.astype(str)).all()
    rows_match = sample_ids.isin(mvals.index.astype(str)).all()
    if not cols_match and rows_match:
        mvals = mvals.T
        print("Transposed methylation matrix to match sample orientation")
    elif not cols_match and not rows_match:
        raise ValueError(
            "Sample IDs do not match methylation matrix dimensions"
        )

    pheno = pheno.set_index(args.sample_id_col)
    mvals = mvals.loc[:, pheno.index]

    design = add_constant(pheno[[args.assoc]].astype(float)).values

    results = []
    for start in range(0, mvals.shape[0], args.chunk_size):
        chunk = mvals.iloc[start : start + args.chunk_size].values.T
        # n_samples x chunk_size
        betas, pvals = run_vectorized(design, chunk)
        cpgs = mvals.index[start : start + args.chunk_size]
        for beta, pval, cpg in zip(betas, pvals, cpgs):
            results.append({"CpG": cpg, "beta": beta, "pvalue": pval})

    out_df = pd.DataFrame(results)
    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, f"ewas_results{args.out_type}")
    out_df.to_csv(out_path, index=False)
    print(f"Results saved to {out_path}")


if __name__ == "__main__":
    main()

