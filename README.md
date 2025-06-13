# EWAS Python Pipeline

This is the Python version of the Epigenome-Wide Association Study (EWAS) pipeline.
It includes scripts for:

- Running EWAS (regression)  
- Stratifying data  
- Annotating CpG sites  
- Creating Manhattan and QQ plots  
- Generating BED files for DMR analysis

## Setup

```bash
pip install -r requirements.txt
```

## CpG Annotation Files

To annotate CpG sites you need annotation tables for the Illumina EPIC array.
You can download `EPIC.hg38.manifest.tsv.gz` and the accompanying SNP key from
the [InfiniumAnnotation project](https://zwdzwd.github.io/InfiniumAnnotation/).
Place the files under a directory named `annotation_files` in the repository
root and rename them as follows:

```bash
mkdir -p annotation_files
curl -L -o annotation_files/EPIC_hg38.tsv.gz \
  https://zwdzwd.github.io/InfiniumAnnotation/files/EPIC.hg38.manifest.tsv.gz
curl -L -o annotation_files/EPIC_snp_key.tsv.gz \
  https://zwdzwd.github.io/InfiniumAnnotation/files/EPIC.hg38.snp.tsv.gz
```

The scripts will look for these filenames inside `annotation_files/` by
default.

## Usage Examples

Assuming your phenotype table (`pheno.csv`) and methylation matrix (`mvals.csv.gz`) are stored in the `data` directory:

### Run EWAS
```bash
python scripts/ewas.py --pheno data/pheno.csv --methyl data/mvals.csv.gz --assoc BMI --out-dir output/ewas
```
Use `--sample-id-col` if your phenotype table uses a different column name for sample IDs. The script will automatically transpose the methylation matrix if samples are found on rows.

### Run EWAS (vectorized)
```bash
python scripts/ewas_fast.py --pheno data/pheno.csv --methyl data/mvals.csv.gz --assoc BMI --out-dir output/ewas
```
This alternative implementation performs the regressions in chunks using NumPy operations to avoid repeatedly copying the design and methylation matrices.

### Stratify the data
```bash
python scripts/stratify.py --pheno data/pheno.csv --methyl data/mvals.csv.gz --stratify sex re --out-dir output/stratified
```
If your phenotype table uses another sample column, specify it with `--sample-id-col`.

### Annotate EWAS results
```bash
python scripts/annotate.py --input-file output/ewas/ewas_results.csv --out-dir output/annotated --assoc BMI --stratified no
```
Use `--anno-file` to provide a different CpG annotation file if needed.
Annotation files such as `EPIC_hg38.tsv.gz` and `EPIC_snp_key.tsv.gz` are expected inside an `annotation_files` directory (see [CpG Annotation Files](#cpg-annotation-files) for download instructions).

### Create a BED file
```bash
python scripts/make_bed.py --results output/annotated/BMI_ewas_annotated_results.csv --assoc BMI --out-dir output/bed --p-threshold 1e-5
```

### Generate Manhattan and QQ plots
```bash
python scripts/plots.py --input-file output/annotated/BMI_ewas_annotated_results.csv --assoc BMI --out-dir output/plots
```
Use `--out-type png` to save the plots as PNG images instead of the default HTML.
The option accepts either `html` or `png`.

Results will be written to the paths supplied with `--out-dir`; the examples above place all output under the `output` directory.
