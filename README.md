# Perturb-Seq workflow

This repository provides a reproducible, Python-driven workflow that starts with
10x Genomics FASTQ files and finishes with cell-level annotations, guide calls,
quality-control summaries, embeddings, and perturbation differential-expression
tables. Python controls every step; Cell Ranger is invoked for the alignment and
UMI-counting step because it understands the 10x feature-barcode chemistry.

> **Assumptions.** The experiment uses 10x gene-expression plus CRISPR Guide
> Capture libraries. The FASTQ sample name (`sample` below) must match the prefix
> produced by `cellranger mkfastq`/`bcl-convert`. Always confirm chemistry, guide
> feature-reference sequences, genome build, and experimental covariates before
> interpreting results.

## 1. Install software

Install Cell Ranger separately and ensure `cellranger` is on `PATH`. Then create
a Python environment (Python 3.10 or newer):

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

For development and tests, also run `python -m pip install -r requirements-dev.txt`.

## 2. Prepare references and FASTQs

Download or build a Cell Ranger transcriptome reference. Create the feature
reference CSV required by Cell Ranger; guide sequences must be the protospacer
sequence as read in the feature-barcode read. For example:

```csv
id,name,read,pattern,sequence,feature_type,target_gene_id,target_gene_name
NTC_1,NTC_1,R2,(BC),GCGAGGTATTCGGCTCCGCG,CRISPR Guide Capture,NTC,NTC
TP53_1,TP53_1,R2,(BC),ACGTACGTACGTACGTACGT,CRISPR Guide Capture,TP53,TP53
```

Keep all lanes for one library together; Cell Ranger discovers lane/read files
from their standard names. Do not concatenate or rename individual FASTQs.

## 3. Configure the run

Copy the documented template and edit paths and thresholds:

```bash
cp config/example.yaml config/run.yaml
python -m perturb_seq validate --config config/run.yaml
```

The `controls` values must exactly match guide `id` values in the feature
reference. Add biological covariates such as batch or donor to the metadata after
counting and before model fitting if the experiment needs them.

## 4. FASTQ to count matrix

Inspect the generated Cell Ranger command first, then execute it:

```bash
python -m perturb_seq count --config config/run.yaml --dry-run
python -m perturb_seq count --config config/run.yaml
```

This writes Cell Ranger output under `work_dir/<run_id>/`. Preserve `web_summary.html`
and `metrics_summary.csv`: they reveal sequencing saturation, mapping, cell calls,
and feature-barcode performance. A failed Cell Ranger command stops the workflow.

If counting was performed elsewhere, set `matrix_h5` in the YAML to its
`filtered_feature_bc_matrix.h5` and skip this step.

## 5. Guide assignment, QC, and exploratory analysis

```bash
python -m perturb_seq analyze --config config/run.yaml
```

The analysis performs these operations in order:

1. Read the filtered gene/feature matrix and retain both gene expression and
   CRISPR Guide Capture counts.
2. Assign a guide only when it has at least `min_guide_umis` and exceeds the
   second guide by `guide_ratio`; otherwise label the cell `unassigned` or
   `multiplet`. This transparent heuristic should be checked against the guide
   UMI distributions for each experiment.
3. Calculate mitochondrial fraction, total UMIs, detected genes, and guide UMIs;
   apply the configured cell and gene QC thresholds.
4. Save raw counts in a layer, library-size normalize, log-transform, choose
   highly variable genes, scale, run PCA/neighbors/UMAP, and optionally Leiden.
5. Compare every sufficiently represented single-guide perturbation with pooled
   non-targeting controls using a cell-level Wilcoxon screen. Report Benjamini–
   Hochberg adjusted p-values and log fold changes. For publication, validate hits
   with replicate-aware pseudobulk or a dedicated Perturb-seq model rather than
   treating cells as independent biological replicates.

## 6. Final outputs

Results are placed in `output_dir`:

| Output | Meaning |
| --- | --- |
| `perturb_seq_processed.h5ad` | Complete annotated analysis object |
| `cell_metadata.csv.gz` | Per-cell QC, guide assignment, cluster, and UMAP values |
| `guide_assignment_summary.csv` | Cell counts per guide/assignment state |
| `qc_summary.json` | Input/output cell numbers and configured thresholds |
| `differential_expression.csv.gz` | Ranked per-guide-versus-control results |
| `figures/qc_violin.png` | Cell QC distributions |
| `figures/umap_guide.png` | UMAP colored by guide assignment |

Archive the edited YAML, software versions (`cellranger --version` and
`python -m pip freeze`), Cell Ranger metrics, and these outputs. Review guide
representation, control behavior, replicate concordance, perturbation target
knockdown, and potential batch effects before drawing biological conclusions.

## Command reference

```bash
python -m perturb_seq --help
python -m perturb_seq validate --config config/run.yaml
python -m perturb_seq count --config config/run.yaml [--dry-run]
python -m perturb_seq analyze --config config/run.yaml
```

## Testing

```bash
pytest
ruff check .
```
