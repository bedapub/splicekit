# Differential splicing (edgeR)

Running edgeR analysis on features (junctions, anchors, exons, genes) is a single command:

```bash
splicekit edgeR              # all feature types
splicekit edgeR junctions    # a single feature type
splicekit edgeR exons
splicekit edgeR anchors
splicekit edgeR genes
```

Donor/acceptor anchor results are then merged back into the corresponding junction results by `splicekit juan` (part of `splicekit process`), so a junction's row also carries its anchors' edgeR statistics.

## Results files

Results are stored in `results/results_edgeR_{feature_type}.tab`, where `feature_type` is one of `genes`, `exons`, `junctions`, `donor_anchors`, `acceptor_anchors`. Only results with `FDR < splicekit.config.edgeR_FDR_thr` are reported (sorted by FDR), each linked to [JBrowse2](jbrowse2.md) via a URL.

To explore all results without the FDR filter, use `results/results_edgeR_{feature_type}_all.tab`.

### General columns

| Column | Example | Description |
|---|---|---|
| `result_id` | `r1` | Integer result identifier, starting at 1. |
| `comparison` | `test_control` | Comparison name, from `annotation/comparisons.tab`. |
| `compound` | `treatment1` | Name of the treatment/compound tested. |
| `feature_id` | `chr1+_17741_17839` | ID of the reported feature: a gene/exon/junction/`[donor,acceptor]_anchor` ID. |
| `chr` | `1` | Chromosome of the feature. |
| `strand` | `+` | Strand of the feature (`+` or `-`). |
| `feature_start` | `17741` | Start of the feature (numerically, start < stop). See [Genomic coordinates](coordinates.md). |
| `feature_stop` | `17839` | Stop of the feature (numerically, stop > start). |
| `feature_length` | `250` | `feature_stop - feature_start + 1`. |
| `gene_id` | `ENSG00000120948` | Ensembl or RefSeq gene ID. |
| `gene_name` | `TARDBP` | Corresponding to `gene_id`. |
| `sum_feature_test` | `1000` | Sum of counts for this feature across all test samples. |
| `sum_feature_control` | `1000` | Sum of counts for this feature across all control samples. |
| `jbrowse_loc` | `3:342321..351243` | Genomic region shown in the JBrowse view. |
| `jbrowse_url` | | Link to the JBrowse view. |
| `logFC` | | Log fold change, from edgeR. |
| `exon.F` | | `exon.F` statistic, from edgeR. |
| `p_value` | | p-value, from edgeR. |
| `fdr` | | False discovery rate, from edgeR. |

### Junction-specific (additional) columns

| Column | Example | Description |
|---|---|---|
| `annotated` | `AA` | Two-letter code `AA`/`AN`/`NA`/`NN`: first letter for the donor site (5' of junction), second for the acceptor site (3' of junction); `A` = touches an annotated exon, `N` = does not. |
| `donor_anchor_id` | | ID of the donor anchor linked to this junction. |
| `acceptor_anchor_id` | | ID of the acceptor anchor linked to this junction. |
| `UTR` | | `first_exon_{start_pos}` if the junction touches any transcript's first exon of the gene. |

### Exon-specific (additional) columns

| Column | Example | Description |
|---|---|---|
| `delta_PSI` | | `test_PSI - control_PSI` (percentage spliced-in). |

Next: [Motif & RNA-binding analysis](motifs.md) runs on the sequences around the regulated features found here.
