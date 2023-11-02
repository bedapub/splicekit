<picture><img src="media/splicekit_logo.png" height="30"/></picture>
## splicekit: an integrative toolkit for splicing analysis from short-read RNA-seq

<b>splicekit</b> is a modular platform for splicing analysis from short-read RNA-seq datasets. The platform also integrates an JBrowse2 instance, [pybio](https://github.com/grexor/pybio) for genomic operations and [scanRBP](https://github.com/grexor/scanRBP) for RNA-protein binding studies. The whole analysis is self-contained (one single folder) and the platform is written in Python, in a modular way.

## Quick start

The easiest way to install splicekit is to simply run:

`$ pip install splicekit`

Note that on some systems, pip is installing the executable scripts under `~/.local/bin`. However this folder is not in the PATH which will result in `command not found` if you try to run `$ splicekit` on the command line. To fix this, please execute:

`export PATH="$PATH:~/.local/bin"`

Another suggestion is to install inside a virtual environment (using `virtualenv`).

<details>
<summary>Installing splicekit directly from the GitHub repository</summary>

```
# move to HOME folder (for example)
cd ~
# clone the repository
git clone git@github.com:bedapub/splicekit.git
# adjust PYTHONPATH and PATH
export PYTHONPATH=$PYTHONPATH:~/splicekit
export PATH=$PATH:~/splicekit/splicekit
```
</details>

<details>
<summary>If you already have aligned reads in BAM files</summary>

All you need is `samples.tab` (note that this is a <b>TAB delimited file</b>) and `splicekit.config` in one folder (check [datasets](datasets) for examples).

You can easily download and prepare the reference genome (e.g. `$ pybio genome homo_sapiens`).

Finally run `$ splicekit process` (inside the folder with `samples.tab` and `splicekit.config`).

Easiest is to check [datasets](datasets) examples to see how the above files look like and also to check scripts if you need to map reads from FASTQ files with `pybio`.
</details>

## Documentation

* [PDF reference manual](https://github.com/bedapub/splicekit/raw/main/docs/splicekit_docs.pdf)
* [Google docs](https://docs.google.com/document/d/15ZRCeK8xyg3klLktZSHZ9k__Xw_BZRn_Q-J4W35JNnQ/edit?usp=sharing) of the above PDF (comment if you like)

## Changelog<a name="changelog"></a>

**v0.5**: released in September 2023
* added rMATS analysis for splicing events
* added Docker container that can be directly imported to singularity via ghcr.io

<details>
<summary><b>v0.4</b>: released in May 2023 (click to show details)</summary>

* added singularity container with all dependencies
* added local integrated JBrowse2
* cluster or desktop runs
* scanRBP and bootstrap analysis of RNA-protein binding
* further development and integration with pybio
* extended documentation of concepts, analysis and results
</details>

<details>
<summary><b>v0.3</b>: released in January 2023 (click to show details)</summary>

* re-coded junction analysis
  * independent junctions parsing from provided bam files
  * master table of all junctions in the samples of the analyzed project, including novel junctions (refseq/ensembl non-annotated)
* clustering by logFC of pairwise-comparisons with dendrogram: junction, exon and gene levels (clusterlogfc module)
* added *first_exon* annotation for junctions touching annotated first exons of transcripts
* extended documentation of concepts, analysis and results
</details>

<details>
<summary><b>v0.2</b>: released in October 2022 (click to show details)</summary>

* software architecture restructure with python modules
* filtering of lowly expressed features by edgeR
* DonJuan analysis (junction-anchor analysis)
* more advanced motif analysis with DREME
* filtering regulated junctions with regulated donors
</details>

<details>
<summary><b>v0.1</b>: released in July 2022 (click to show details)</summary>

* initial version of splicekit
* parsing of junction and exon counts
* computing edgeR analysis from count tables and producing a results file with direct links to JBrowse2
* basic motif analysis
</details>
