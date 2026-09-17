# Additional analyses

Beyond the core edgeR / motif / juDGE pipeline, splicekit ships several further analyses, each runnable as its own `splicekit` command and included in `splicekit process`.

## Promiscuity: `splicekit promisc`

```bash
splicekit promisc
```

Promiscuity analysis looks at genes and how many of their junctions change across conditions — genes with many independently-regulated junctions ("promiscuous" splicing) versus genes regulated through a single dominant junction.

## Cluster of pairwise logFC: `splicekit clusterlogfc`

```bash
splicekit clusterlogfc
```

Clusters comparisons by the logFC of their significant (FDR < 0.05) changes, computed separately at the junction, exon and gene level. This groups comparisons/treatments that produce similar splicing/expression signatures.

## JUNE: junction-event analysis

```bash
splicekit june
```

JUNE (junction-events) classifies pairs/sets of overlapping junctions into splicing event types — e.g. skipped exons and mutually exclusive exons — by comparing their shared and differing donor/acceptor coordinates, rather than looking at single junctions or exons in isolation.

## rMATS

```bash
splicekit rmats
```

Runs [rMATS](http://rnaseq-mats.sourceforge.net/) on the project's BAM files as an independent, complementary splicing-event caller (skipped exons, alternative 5'/3' splice sites, mutually exclusive exons, retained introns) alongside splicekit's own junction/exon-based analysis.

## DEXSeq

```bash
splicekit dexseq              # all feature types
splicekit dexseq junctions    # a single feature type
```

An alternative to `splicekit edgeR` for differential feature usage, using [DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html) instead of edgeR. Requires `dexseq_scripts` to be set in `splicekit.config` to the path of the DEXSeq R scripts.

## DonJuAn: `splicekit juan`

```bash
splicekit juan
```

Merges donor/acceptor anchor edgeR results back into the junction results (`results/results_edgeR_junctions.tab`), so each junction's row also reports its anchors' differential usage statistics. Runs automatically as part of `splicekit process`, right after `splicekit edgeR`.

## HTML report: `splicekit report`

```bash
splicekit report
```

Assembles the project's results (edgeR, JUNE, and more) into a browsable HTML report under `report/`, served together with [JBrowse2](jbrowse2.md) when you run `splicekit web`.
