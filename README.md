<picture><img src="media/splicekit_logo.png" height="30"/></picture>
## splicekit: an integrative toolkit for splicing analysis from short-read RNA-seq

<b>splicekit</b> is a modular platform for splicing analysis from short-read RNA-seq datasets. The platform also integrates an JBrowse2 instance, [pybio](https://github.com/grexor/pybio) for genomic operations and [scanRBP](https://github.com/grexor/scanRBP) for RNA-protein binding studies. The whole analysis is self-contained (one single folder) and the platform is written in Python, in a modular way.

Check a short video presentation about splicekit (poster) at ECCB 2023 on Youtube:

[<img src="media/splicekit_youtube.jpg" width=300>](https://youtu.be/P1m73usZ3lc?si=HBJxWOkUajObFpu1)

## Quick start

Since version 0.7, splicekit is a **Snakemake** pipeline, and there is also a Conda environment yaml file.

```
git clone git@github.com:bedapub/splicekit.git    # clone rep
cd splicekit                                      # change working directory

micromamba -y create -f splicekit.yaml            # create conda env
micromamba activate splicekit                     # activate env
./install.sh                                      # install dependencies
pip install .                                     # install splicekit

cd datasets/GSE126543                             # move to sample folder
./1_download.sh                                   # download sample FASTQs
pybio homo_sapiens                                # human genome

./run_snakemake_local.sh --configfile config.yaml # run snakemake (local)
./run_snakemake_slurm.sh --configfile config.yaml # OR run snakemake SLURM
```

After snakemake finishes, you can explore results interactively running `splicekit web` and follow instructions on how to open the html reports in your browser.

<details>
<summary>Installing splicekit directly from the GitHub repository</summary>

```
pip install git+https://github.com/bedapub/splicekit.git@main
```
</details>

<details>
<summary>If you already have aligned reads in BAM files</summary>

All you need is `samples.tab` (note that this is a <b>TAB delimited file</b>) and `splicekit.config` in one folder (check [datasets](datasets) for examples).

You can easily download and prepare the reference genome (e.g. `$ pybio genome homo_sapiens`).

Finally run `./run_snakemake_[local/slurm].sh --configfile config.yaml` (inside the folder with `samples.tab` and `splicekit.config`).

Easiest is to check [datasets](datasets) examples to see how the above files look like and also to check scripts if you need to map reads from FASTQ files with `pybio`.
</details>

## Documentation

* [PDF reference manual](https://github.com/bedapub/splicekit/raw/main/docs/splicekit_docs.pdf)
* [Google docs](https://docs.google.com/document/d/15ZRCeK8xyg3klLktZSHZ9k__Xw_BZRn_Q-J4W35JNnQ/edit?usp=sharing) of the above PDF (comment if you like)

## Changelog<a name="changelog"></a>

**v0.7**: released in February 2025

* [Snakemake](Snakefile) version
* Conda [splicekit.yaml](splicekit.yaml) for environment setup

<details>
<summary>Past change notes (click to view)</summary>
<br>
<b>v0.6</b>: released in April 2024

* updated reports
* JUNE analysis (junction-events to classify skipped and mutually exclusive exons)

<b>v0.4.9</b>: released in November 2023

* added rMATS analysis for splicing events
* added Docker container that can be directly imported to singularity via ghcr.io
* fixed dependencies
* other small fixes

<b>v0.4</b>: released in May 2023

* added singularity container with all dependencies
* added local integrated JBrowse2
* cluster or desktop runs
* scanRBP and bootstrap analysis of RNA-protein binding
* further development and integration with pybio
* extended documentation of concepts, analysis and results

<b>v0.3</b>: released in January 2023 (click to show details)

* re-coded junction analysis
  * independent junctions parsing from provided bam files
  * master table of all junctions in the samples of the analyzed project, including novel junctions (refseq/ensembl non-annotated)
* clustering by logFC of pairwise-comparisons with dendrogram: junction, exon and gene levels (clusterlogfc module)
* added *first_exon* annotation for junctions touching annotated first exons of transcripts
* extended documentation of concepts, analysis and results

<b>v0.2</b>: released in October 2022

* software architecture restructure with python modules
* filtering of lowly expressed features by edgeR
* DonJuan analysis (junction-anchor analysis)
* more advanced motif analysis with DREME
* filtering regulated junctions with regulated donors

<b>v0.1</b>: released in July 2022

* initial version of splicekit
* parsing of junction and exon counts
* computing edgeR analysis from count tables and producing a results file with direct links to JBrowse2
* basic motif analysis

</details>

## Citing and Contact<a name="citation"></a>

If you find **splicekit** useful in your work and research, please cite:

Rot, G., Wehling, A., Schmucki, R., Berntenis, N., Zhang, J. D., & Ebeling, M. (2024)<br>
[splicekit : an integrative toolkit for splicing analysis from short-read RNA-seq](https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae121/7735317)<br>
Bioinformatics Advances, 4(1). https://doi.org/10.1093/bioadv/vbae121

In case of questions, issues and other ideas, please use the <a href='https://github.com/bedapub/splicekit/issues'>GitHub Issues</a> or write directly to <a href='mailto:gregor.rot@gmail.com'>Gregor Rot</a>.
