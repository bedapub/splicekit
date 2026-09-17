# Quick Start

To run splicekit you need:

1. **A reference genome**, downloaded and processed with `pybio` (installed automatically as a splicekit dependency):
    ```bash
    pybio genome homo_sapiens        # or: pybio search species, for other species
    ```
2. **Aligned reads in BAM format**, one file per sample. You can align FASTQ files yourself with STAR, or reuse a dataset's mapping script, e.g. [datasets/GSE221868/2_map.sh](https://github.com/bedapub/splicekit/blob/main/datasets/GSE221868/2_map.sh), which downloads the reference genome with pybio and aligns with STAR.
3. **`samples.tab`** — one line per sample, TAB delimited, connecting each `sample_id` to its `treatment_id`. See [Sample annotation](samples.md) and the [example samples.tab](https://github.com/bedapub/splicekit/blob/main/datasets/GSE126543/samples.tab).
4. **`splicekit.config`** — reference genome, BAM folder and the other core parameters. See [Configuration](configuration.md) and the [example splicekit.config](https://github.com/bedapub/splicekit/blob/main/datasets/GSE126543/splicekit.config).
5. **`config.yaml`** — per-rule Snakemake resources (cores/memory/time). Copy the [template config.yaml](https://github.com/bedapub/splicekit/blob/main/config.yaml) into your project folder and adjust it to your cluster/machine.

The [datasets](https://github.com/bedapub/splicekit/tree/main/datasets) folder has four complete examples, each with its own scripts to download and process a public RNA-seq dataset from scratch.

## Running the pipeline

With `samples.tab`, `splicekit.config` and `config.yaml` in your project folder, run the whole pipeline with Snakemake:

```bash
cd datasets/GSE126543                               # example project folder
./1_download.sh                                     # download sample FASTQs
pybio homo_sapiens                                   # reference genome

./run_snakemake_local.sh --configfile config.yaml    # run locally
# or:
./run_snakemake_slurm.sh --configfile config.yaml    # submit jobs to SLURM
```

`run_snakemake_slurm.sh` submits each Snakemake rule as its own SLURM job (via `snakemake-executor-plugin-cluster-generic`), sized per-rule from `config.yaml`.

Once it finishes, explore the results:

```bash
splicekit web
```

This starts a single local web server serving both the HTML report (`http://<host>:8007/report`) and the [JBrowse2](jbrowse2.md) genome browser.

!!! note
    If you already have BAM files and want to skip Snakemake, you can run splicekit directly with `splicekit process` inside a folder containing `samples.tab` and `splicekit.config` — this runs the same analysis steps sequentially on a single machine. See [Command-line reference](cli.md).

## Next steps

- [Configuration](configuration.md) — every `splicekit.config` and `config.yaml` parameter.
- [Sample annotation](samples.md) — how `samples.tab` becomes `annotation/comparisons.tab`.
- [Features & count tables](features.md) — the four feature types and their count files.
- [Differential splicing (edgeR)](edgeR.md) — running and reading edgeR results.
