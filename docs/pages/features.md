# Features & count tables

Running `splicekit features` creates count tables for junctions, anchors, exons and genes from the project's BAM files.

## What are features?

splicekit operates on 4 types of features: **junctions, anchors, exons and genes**. All feature IDs share the same format: `chrstrand_start_stop`, e.g. `chrX-_154371360_154374505`. See [Genomic coordinates](coordinates.md) for how coordinates are reported across splicekit.

Junctions are detected directly from BAM files (independent of any pre-existing gene model), and reported in `reference/junctions.tab` together with a *donor anchor* and *acceptor anchor* — by default the 15nt regions flanking the junction's start and stop. These anchor regions are turned into `reference/donor_anchors.gtf` / `reference/acceptor_anchors.gtf` and quantified with **featureCounts**, alongside exon- and gene-level counts. See [File formats](fileformats.md) for the exact column layouts.

## Feature data files

Each individual sample gets one file (table) per feature type, under `data/sample_{feature_type}_data/`, listing every feature and its count in that sample.

Example: `data/sample_exons_data/sample_99.tab`

```
GeneID  Start     End       Length  Symbol  1_test  2_test  3_control  4_control
1       58347029  58347353  325     A1BG    42      31      109        75
1       58347640  58350370  2731    A1BG    0       0       3          1
1       58350651  58351391  741     A1BG    0       0       10         1
```

Next step: [Differential splicing (edgeR)](edgeR.md), which turns these per-sample count tables into per-comparison differential usage results.
