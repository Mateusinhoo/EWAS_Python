"""Annotate EWAS results with genomic information."""

import argparse
import os
import pandas as pd

def main() -> None:
    """Entry point for command line usage."""

    parser = argparse.ArgumentParser(
        description="Add annotation data to EWAS results"
    )
    parser.add_argument(
        "--input-file",
        "-i",
        required=True,
        help="EWAS result file (CSV or TSV)",
    )
    parser.add_argument(
        "--out-dir", required=True, help="Directory to save annotated results"
    )
    parser.add_argument(
        "--stratified",
        choices=["yes", "no"],
        default="no",
        help="Was the analysis stratified",
    )
    parser.add_argument("--assoc", required=True, help="Association variable name")
    parser.add_argument("--out-type", default=".csv", help="Output file extension")
    parser.add_argument(
        "--anno-file",
        default="annotation_files/EPIC_hg38.tsv.gz",
        help="Path to annotation file",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input_file, sep=None, engine="python")

    anno_cpg = pd.read_csv(args.anno_file, sep="\t")

    annotated = df.merge(anno_cpg, on="CpG", how="left")

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(
        args.out_dir, f"{args.assoc}_ewas_annotated_results{args.out_type}"
    )
    annotated.to_csv(out_path, index=False)
    print(f"Annotated results saved to {out_path}")


if __name__ == "__main__":
    main()
