# splicekit: splicing analysis from short-read RNA-seq

![splicekit](assets/splicekit_logo.png)

**splicekit** is a modular, integrative platform for splicing analysis of short-read RNA-seq data. Starting from aligned reads (BAM files) and a sample annotation table, it defines test-vs-control comparisons, builds per-feature count tables (junctions, anchors, exons, genes), and runs a battery of splicing analyses — differential feature usage with edgeR, motif and RNA-protein binding enrichment, junction-vs-gene expression comparisons and more — all self-contained in a single project folder. It integrates [pybio](https://github.com/grexor/pybio) for genome operations, [scanRBP](https://github.com/grexor/scanRBP) for RNA-protein binding, and ships its own [JBrowse2](https://jbrowse.org/jb2/) instance for browsing results.

```bash
# create and activate the conda/micromamba environment, install splicekit
micromamba create -y -f splicekit.yaml
micromamba activate splicekit
./install.sh
pip install .

# run the full pipeline with Snakemake
./run_snakemake_local.sh --configfile config.yaml
```

## What's included

- **Comparisons from a sample sheet** — define test/control comparisons straight from `samples.tab`, optionally grouped or separated by extra columns. See [Sample annotation](samples.md).
- **Feature count tables** — junctions, anchors, exons and genes, built from BAM files. See [Features & count tables](features.md).
- **Differential splicing with edgeR** — per-feature differential usage results, linked directly to JBrowse2. See [Differential splicing (edgeR)](edgeR.md).
- **Motif & RNA-protein binding analysis** — donor/acceptor motif logos, DREME enrichment and scanRBP binding analysis. See [Motif & RNA-binding analysis](motifs.md).
- **juDGE plots** — junction logFC vs. gene logFC, distinguishing splicing modifiers from expression modifiers. See [juDGE plots](judge.md).
- **Promiscuity, clustering, JUNE and rMATS analyses** — see [Additional analyses](additional-analyses.md).
- **Integrated JBrowse2 + HTML report** — one local web server for both. See [Exploring results](jbrowse2.md).

## Where to start

New to splicekit? Read [Installation](installation.md) and then [Quick Start](quickstart.md) — together they take you from a fresh checkout to a running pipeline on example data in a few commands. Everything else in these docs is reference material for the individual analysis steps, configuration parameters and file formats.
