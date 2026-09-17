# Dependencies

## Conda/micromamba environment (`splicekit.yaml`)

Installed by `micromamba create -y -f splicekit.yaml`:

[pigz](https://zlib.net/pigz/), [deeptools](https://deeptools.readthedocs.io/), [samtools](http://www.htslib.org/), [Snakemake](https://snakemake.readthedocs.io/), R (`r-base`, `r-locfit`), [STAR](https://github.com/alexdobin/STAR), [rMATS](http://rnaseq-mats.sourceforge.net/), [Node.js](https://nodejs.org/), Ghostscript, [subread](http://subread.sourceforge.net/) (`featureCounts`, > 2.0.6), [MEME](https://meme-suite.org/meme/) (> 5.5.1, for DREME), Perl's `cpanminus`, and the `snakemake-executor-plugin-cluster-generic` pip package (used by `run_snakemake_slurm.sh` for SLURM submission).

## Installed by `install.sh`

- R packages: `BiocManager`, `data.table`, `statmod`, `R.utils`, and Bioconductor's [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html).
- Perl modules for the JBrowse2/SOAP tooling (`XML::Compile::*`, `Log::Log4perl`, `Math::CDF`, `JSON`, and others).
- [`@jbrowse/cli`](https://www.npmjs.com/package/@jbrowse/cli) via `npm install -g`, used to set up the local [JBrowse2](jbrowse2.md) instance.

## Python dependencies (installed by `pip install .`)

[Levenshtein](https://pypi.org/project/Levenshtein/), [logomaker](https://logomaker.readthedocs.io/), [plotly](https://plotly.com/python/), `python-dateutil`, [pybio](https://grexor.github.io/pybio/), [scanRBP](https://grexor.github.io/scanRBP/), [pysam](https://pysam.readthedocs.io/), [numpy](https://numpy.org/), `psutil`, `beautifulsoup4`, `requests` and `rangehttpserver`.

## Optional

[Singularity](https://sylabs.io/singularity/) — only needed if you set `container = "singularity run docker://ghcr.io/bedapub/splicekit:main"` in `splicekit.config` instead of installing the conda environment directly. See [Installation: Container](installation.md#container).
