# File formats

## `reference/junctions.tab`

Contains all junctions detected across every sample in the project. Only junctions that could be annotated to a gene are reported — including "novel" junctions that don't touch a RefSeq/Ensembl-annotated exon, as long as the junction's start and stop fall inside an annotated gene (see the `annotated` column).

| Column | Example | Description |
|---|---|---|
| `junction_id` | `chr1+_17741_17839` | Unique ID: `chrstrand_start_stop`. |
| `donor_anchor_id` | `chr1+_17725_17740` | Matching donor anchor ID — by default the 15nt region upstream of the junction start. |
| `acceptor_anchor_id` | `chr1+_17840_17855` | Matching acceptor anchor ID — by default the 15nt region downstream of the junction stop. |
| `gene_id` | `ENSG00000120948` | Ensembl or RefSeq gene ID. A junction can be non-annotated (`annotated != "AA"`) but still assigned to a gene, meaning its start/stop fall inside the gene. |
| `gene_name` | `TARDBP` | Corresponding to `gene_id`. |
| `chr` | `1` | Chromosome. |
| `strand` | `+` | `+` or `-`. |
| `annotated` | `AA` | Two-letter code `AA`/`AN`/`NA`/`NN` — see [Genomic coordinates](coordinates.md) and [edgeR results](edgeR.md#junction-specific-additional-columns). |
| `count` | `553` | Raw read count across all samples in the project supporting this junction. |

## `reference/donor_anchors.gtf` and `reference/acceptor_anchors.gtf`

GTF files generated from all donor/acceptor anchors in `reference/junctions.tab`. Used by **featureCounts** to build anchor count tables across the project's samples.

## `results/results_edgeR_{feature_type}.tab`

See [Differential splicing (edgeR): Results files](edgeR.md#results-files) for the full column reference (general columns, plus junction- and exon-specific additions).

## `data/sample_{feature_type}_data/*.tab`

See [Features & count tables: Feature data files](features.md#feature-data-files).

## `annotation/comparisons.tab`

See [Sample annotation: Comparisons](samples.md#comparisons).
