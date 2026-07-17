# Perturb-Seq Workflow

Auditable preprocessing for pooled Perturb-seq experiments with multiple treatment
conditions and expected same-target dual-guide designs. The workflow creates a true
unfiltered ingestion snapshot, explicit QC decisions, filtered count and normalized
layers, synchronized guide counts, and sparse downstream exports.

## Input contract

Copy `config.example.yaml` and provide at least one input per condition. Each
preprocessing input requires `input_id`, `condition`, and a Cell Ranger
`matrix_h5`. Barcodes become `<input_id>:<original_barcode>`, so identical 10x
barcodes in different libraries cannot collide. FASTQ and reference fields are not
required for matrix-only preprocessing.

The guide annotation is a TSV/CSV with these columns:

| column | meaning |
|---|---|
| `guide_id` | feature ID exactly matching the Cell Ranger CRISPR feature |
| `target` | target gene, or a control label |
| `guide_type` | targeting, non-targeting, etc. |
| `control_status` | `non_targeting` for NTCs; otherwise `targeting` |
| `expected_pair` | partner guide ID for an intended dual-guide pair; blank otherwise |

## Run

```bash
python -m pip install -e .
perturb-seq preprocess config.yaml
```

The pipeline intentionally does not run differential expression. When biological
replicates are absent, drug, knockdown, and interaction results must be described as
exploratory rather than treated as replicated inference.

## Output contract

Outputs are stage-specific:

* `00_ingestion/`: unfiltered RNA and guide H5ADs plus barcode/feature alignment report.
* `01_qc/`: annotated pre-QC object, QC-filtered object, cell decision log, and gene log.
* `02_final/final_expression.h5ad`: `X` is log1p normalized; `layers["counts"]`
  contains filtered integer counts; `layers["normalized"]` contains pre-log values.
* `02_final/guide_counts.h5ad`: guide counts in exactly the same retained-cell order.
* `02_final/`: gzipped Matrix Market counts/log1p matrices, barcodes, gene and cell
  metadata, guide assignments, and a checksum-bearing `matrix_manifest.json`.

The assignment classes are `no_guide_or_unassigned`, `non_targeting_control`,
`single_guide_targeting`, `expected_same_target_dual`, `ambiguous_assignment`,
`unexpected_multi_target`, and `guide_multiplet`. Intended pairs are determined from
the guide map, not from an inappropriate top/second-guide dominance ratio.

## QC and interpretation

Cell metrics include RNA UMIs, detected genes, mitochondrial and ribosomal fractions,
guide burden/multiplicity, and explicit doublet fields. Every threshold produces a
pass/fail column before subsetting. The default doublet fields say `not_run`; install
and integrate a validated scorer for production rather than silently claiming that
doublet detection occurred. Gene metadata records prevalence, group support, target
status, test eligibility, annotation eligibility, and the filter reason.

Run tests with `pytest`.
