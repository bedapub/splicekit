# Installation

Since v0.7, **splicekit** is a [Snakemake](https://snakemake.readthedocs.io/) pipeline with a Conda/micromamba environment.

```bash
# clone the repository
git clone git@github.com:bedapub/splicekit.git
cd splicekit

# create and activate the conda environment
micromamba create -y -f splicekit.yaml
micromamba activate splicekit

# install remaining (non-conda) dependencies: R/edgeR, Perl modules, jbrowse-cli
./install.sh

# install splicekit itself
pip install .
```

!!! note
    `install.sh` installs R packages (`edgeR` via BiocManager), Perl SOAP/XML modules, and `@jbrowse/cli` via npm — these aren't packaged as conda dependencies, so run it once per environment.

The `splicekit.yaml` environment brings in Snakemake, STAR, samtools, subread (`featureCounts`), MEME (for DREME), rMATS and the `snakemake-executor-plugin-cluster-generic` plugin used for SLURM submission. See [Dependencies](dependencies.md) for the full list.

## Installing just the Python package

If you already have the environment's tools on your `PATH` (or are only using splicekit's Python API / running individual `splicekit` CLI steps by hand), you can install the package on its own:

```bash
pip install splicekit
```

or, from this repository directly:

```bash
pip install git+https://github.com/bedapub/splicekit.git@main
```

!!! note
    On some systems, **pip** installs the executable scripts under `~/.local/bin`. If this folder is not in your `PATH`, running `splicekit` will fail with `command not found`. Fix this with `export PATH="$PATH:~/.local/bin"` (add it to your `~/.profile` to persist across logins). Another option is to install inside a virtual environment (using [virtualenv](https://virtualenv.pypa.io/en/latest/)).

## Container

If you'd rather not install dependencies directly on the machine or cluster, set `container = "singularity run docker://ghcr.io/bedapub/splicekit:main"` in `splicekit.config`. splicekit will then run its non-Python steps through that imported Docker image (pybio and scanRBP are already installed as regular pip dependencies regardless of the container setting). See [Configuration](configuration.md#processing-parameters).
