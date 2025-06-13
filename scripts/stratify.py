"""Split phenotype and methylation data into stratified groups."""

import argparse
import os

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Stratify phenotype and methylation data")
    parser.add_argument("--pheno", required=True, help="Path to phenotype data")
    parser.add_argument("--methyl", required=True, help="Path to methylation matrix")
    parser.add_argument("--stratify", nargs="+", required=True, help="List of variables to stratify by")
    parser.add_argument("--out-dir", default="results/stratified")
    parser.add_argument(
        "--sample-id-col",
        default="sampleID",
        help="Column name for sample identifiers in the phenotype file",
    )
    args = parser.parse_args()

    pheno = pd.read_csv(args.pheno)
    mvals = pd.read_csv(args.methyl, index_col=0)

    if args.sample_id_col not in pheno.columns:
        raise ValueError(f"Phenotype file must contain a '{args.sample_id_col}' column")

    pheno = pheno.set_index(args.sample_id_col)
    mvals = mvals.loc[:, pheno.index]

    groups = pheno.groupby(args.stratify)

    os.makedirs(args.out_dir, exist_ok=True)

    for group_vals, sub_pheno in groups:
        group_name = (
            "_".join(str(v) for v in group_vals) if isinstance(group_vals, tuple) else str(group_vals)
        )
        sub_mvals = mvals.loc[:, sub_pheno.index]

        pheno_path = os.path.join(args.out_dir, f"{group_name}_pheno.csv")
        mvals_path = os.path.join(args.out_dir, f"{group_name}_mvals.csv")

        sub_pheno.to_csv(pheno_path)
        sub_mvals.to_csv(mvals_path)
        print(f"Saved: {group_name}")


if __name__ == "__main__":
    main()
