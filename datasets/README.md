# splicekit datasets

This folder contains complete splicekit dataset/run examples. In each of the dataset subfolders, you can find simple scripts for:

* downloading FASTQ files of public datasets (1_download.sh)
* mapping reads to the reference genome (2_map.sh)
* running splicekit analysis (3_splicekit.sh)

The datasets are:

* [GSE126543](GSE126543)
* [GSE182150](GSE182150)
* [GSE221868](GSE221868)
* [PRJEB42763](PRJEB42763)

## How to run splicekit on the datasets?

Clone this repository, go into one of the dataset folders and:

* run `1_download.sh` script which will download all the required fastq files
* run `2_map.sh` script which will map the reads to the reference genome with pybio
* run `3_splicekit.sh` which will simply do `splicekit process` and run all splicekit analysis on the data in the folder

## What are depencendies to run examples?

Dependencies needed to run the dataset scripts:
* Install [splicekit](https://github.com/bedapub/splicekit)
  * easiest is to simply run: `pip install splicekit`
* downloading FASTQ files from the archives: fastq-dump
  * [Install SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)
* mapping sequence reads to the reference genome: STAR
  * [Install STAR](https://github.com/alexdobin/STAR)
* downloading Ensembl latest human genome assembly and annotation: pybio
  * Install pybio with `pip install pybio`; download and prepare genome for mapping with STAR: `pybio genome homo_sapiens`
* running splicekit on the resulting bam files: splicekit
  * Simply put the dataset `splicekit.config` and `samples.tab` files in an empty folder and run `splicekit process` inside the folder