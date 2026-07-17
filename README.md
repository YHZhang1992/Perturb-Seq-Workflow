# Perturb-Seq Workflow

Auditable, five-stage analysis for pooled Perturb-seq experiments with multiple
treatment conditions and expected same-target dual-guide designs. Each stage has a
small Python API so it can be run separately, inspected, and resumed.

## Analysis stages

1. **Preprocessing** (`perturb_seq.preprocessing`): optionally run Cell Ranger
   from FASTQs with `count_fastqs`, then turn 10x feature matrices into the QC'd
   expression matrix with `prepare_expression_matrix`.
2. **Metadata clean-up** (`perturb_seq.metadata`): standardize spelling and case
   variants with `clean_metadata` and annotate recognized vehicle/control versus
   treatment conditions.
3. **Pseudobulk and association** (`perturb_seq.pseudobulk`): sum raw counts by
   biological sample, condition, and perturbation with `make_pseudobulk`, then run
   an explicitly named treatment-versus-reference comparison with
   `association_test`.
4. **Functional analysis** (`perturb_seq.functional`): use `overrepresentation`
   to test significant association genes against user-provided gene sets.
5. **Visualization** (`perturb_seq.visualization`): create a reproducible volcano
   plot with `plot_volcano`. Install the `analysis` extra for plotting support.

`perturb_seq.analyses` re-exports all of these public functions for convenient
notebook imports; implementation remains separated by stage.

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
python -m pip install -e '.[analysis]'
perturb-seq preprocess config.yaml
```

For FASTQ inputs, run Cell Ranger before matrix preprocessing:

```python
from perturb_seq.analyses import count_fastqs

matrix_h5 = count_fastqs(
    {
        "input_id": "rapamycin_rep1",
        "fastqs": "data/fastqs/rapamycin_rep1",
        "transcriptome": "references/refdata-gex-GRCh38",
        "sample": "rapamycin_rep1",
        "localcores": 8,
    },
    "cellranger_runs",
)
```

Put the returned path in that input's `matrix_h5` field. Cell Ranger must be
installed separately and available on `PATH`.

The remaining stages can then consume the final H5AD without rerunning Cell
Ranger or QC:

```python
import anndata as ad
from perturb_seq.analyses import (
    association_test,
    clean_metadata,
    make_pseudobulk,
    overrepresentation,
    plot_volcano,
)

cells = ad.read_h5ad("results/02_final/final_expression.h5ad")
cells = clean_metadata(cells, {"rapa": "Rapamycin", "dmso": "DMSO"})

# input_id should identify a biological replicate, not merely a treatment label.
pseudobulk = make_pseudobulk(cells)
associations = association_test(
    pseudobulk, condition="Rapamycin", reference="DMSO"
)
associations.to_csv("results/03_association/rapamycin_vs_dmso.tsv", sep="\t", index=False)

gene_sets = {"example_pathway": ["GENE1", "GENE2", "GENE3"]}
enrichment = overrepresentation(associations, gene_sets)
enrichment.to_csv("results/04_functional/enrichment.tsv", sep="\t", index=False)
plot_volcano(associations, "results/05_visualization/volcano.png")
```

Association testing uses Welch's test on log-CPM values when both conditions
have at least two pseudobulk samples and reports Benjamini-Hochberg adjusted
p-values. With fewer replicates, effect sizes are still returned but p-values
are set to 1; these outputs are descriptive rather than inferential. Provide
one distinct `input_id` per true biological replicate and avoid treating cells
as independent replicates.

When biological replicates are absent, drug, knockdown, and interaction results
must be described as exploratory rather than treated as replicated inference.

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
