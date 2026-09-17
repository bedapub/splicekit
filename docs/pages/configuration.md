# Configuration

splicekit reads parameters from two files in your project folder:

- **`splicekit.config`** — a Python file (each line is `exec`'d) describing the study, sample annotation, genome and analysis parameters. Start from the [template](https://github.com/bedapub/splicekit/blob/main/splicekit/splicekit.config.template) or one of the [datasets](https://github.com/bedapub/splicekit/tree/main/datasets) examples.
- **`config.yaml`** — a Snakemake config file describing per-rule compute resources (cores/memory/time). Start from the [template config.yaml](https://github.com/bedapub/splicekit/blob/main/config.yaml).

## splicekit.config

### Core parameters

**`study_name`**, default = `"descriptive short study name"`
Descriptive short study name, arbitrary string describing the study.

**`library_type`**, default = `"paired-end"`
Either `"paired-end"` or `"single-end"`.

**`library_strand`**, default = `"NONE"`
Possible values:

- `"SECOND_READ_TRANSCRIPTION_STRAND"`
- `"FIRST_READ_TRANSCRIPTION_STRAND"`
- `"SINGLE_STRAND"`
- `"SINGLE_REVERSE"`
- `"NONE"`

For unstranded data, use `"NONE"`. For paired-end stranded data, the most common value is `"SECOND_READ_TRANSCRIPTION_STRAND"`, meaning the second read of the pair maps in the transcript direction and the first read maps in the reverse direction. For stranded single-end sequencing, `"SINGLE_STRAND"` means the reads map in the transcript direction, and `"SINGLE_REVERSE"` means the reads map on the opposite strand of the transcripts.

**`edgeR_FDR_thr`**, default = `0.05`
FDR threshold used to filter `results/results_edgeR_{feature_type}.tab`. See [Differential splicing (edgeR)](edgeR.md).

**`dexseq_scripts`**, default = `""`
Path to DEXSeq R scripts, if you want to run `splicekit dexseq` as an alternative to `splicekit edgeR`.

### Sample annotation parameters

**`sample_column`**, default = `"sample_id"`
Which column in `samples.tab` holds the sample IDs. BAM files are then expected to be named `sample_id.bam` (see `bam_path`/`bam_column` below).

**`treatment_column`**, default = `"treatment_id"`
Which column in `samples.tab` defines treatment (test and control labels).

**`control_name`**, default = `"DMSO"`
The text identifying control samples in the treatment column. Non-control samples are compared against these.

**`separate_column`**, default = `""`
Group samples by this column and only compare within a group (empty = don't separate). Useful when, say, samples come from two tissues and you only want within-tissue comparisons: add a tissue column to `samples.tab` and set `separate_column` to its name.

**`group_column`**, default = `""`
Only compare a test sample against controls from the same "domain" (e.g. the same sequencing plate).

### Genome parameters

**`species`**, default = `None`
Genome species, resolved through `pybio`/Ensembl. Any species pybio knows about, or a custom species you registered with pybio via a FASTA+GTF pair.

**`genome_version`**, default = `None`
Ensembl genome version, or your custom genome version. `None` takes the latest Ensembl release.

**`bam_path`**, default = `""`
Folder where BAM files are stored, expected to contain `{sample_id}.bam` for every sample. Absolute, or relative to the project folder.

**`bam_column`**, default = `"bam_file"`
Alternative to `bam_path`: the name of a column in `samples.tab` holding the full BAM path for each sample individually (useful when BAMs live in per-sample subfolders rather than one flat folder). At least one of `bam_path` or `bam_column` must resolve to something — splicekit exits with an error otherwise.

### scanRBP parameters

**`scanRBP`**, default = `True`
Should splicekit run the scanRBP RNA-protein binding analysis as part of `splicekit motifs`? See [Motif & RNA-binding analysis](motifs.md).

**`protein`**, default = `"K562.TARDBP.0"`
ID of the protein PWM to scan regulated sequences against. List available IDs with `scanRBP search <term>`, e.g.:

```
$ scanRBP search hnRNPA
scan_id           protein  tissue description
HNRNPA1.HepG2.00  HNRNPA1  HepG2  heterogeneous nuclear ribonucleoprotein A1
HNRNPA1.HepG2.01  HNRNPA1  HepG2  heterogeneous nuclear ribonucleoprotein A1
...
HNRNPA1.K562.00   HNRNPA1  K562   heterogeneous nuclear ribonucleoprotein A1
```

**`protein_label`**, default = `"tdp43"`
Short descriptive label for the protein, used in file names and plot titles.

### Processing parameters

**`platform`**, default = `"desktop"`
`"desktop"`, `"cluster"` (LSF/`bsub`) or `"SLURM"`. When running under Snakemake, job submission is instead handled by `run_snakemake_local.sh`/`run_snakemake_slurm.sh` and `config.yaml`.

**`container`**, default = `""`
Empty string assumes all software dependencies are installed on the machine/cluster. Set to `"singularity run docker://ghcr.io/bedapub/splicekit:main"` to run non-Python steps from the provided container image instead (pybio and scanRBP are installed via pip regardless).

**`edgeR_memory`**, default = `"16GB"`
Memory reserved for edgeR cluster jobs (only relevant on `"cluster"`/`"SLURM"` platforms; under Snakemake, use `config.yaml` instead).

### Visualization and labeling parameters

**`short_names`**, default = `[]`
By default, no names are shortened or replaced in results or plots. Provide a list of triples to replace/shorten long strings:

```python
# example splicekit.config short_names parameter
# replace cell_line_A with A, cell_line_B with B (only on an exact/"complete" match)
short_names = [("cell_line_A", "A", "complete"), ("cell_line_B", "B", "complete")]
```

Use `"partial"` instead of `"complete"` to also replace the string when it occurs inside a larger one (e.g. `...cell_line_A...` → `...A...`).

## config.yaml

`config.yaml` sizes the Snakemake rules — how many cores, how much memory and how much walltime each step gets, whether running locally or submitting to SLURM.

```yaml
mapping:
  perform_mapping: True  # should splicekit map FASTQ -> BAM with pybio/STAR? (True/False)
  alignIntronMax: 0      # STAR max intron size; 0 lets STAR derive it from its default window parameters

defaults:
  cores: 1
  mem: 4            # GB
  time: "01:00:00"

map_fastq_single:   { cores: 8, mem: 8, time: "02:00:00" }
map_fastq_paired:   { cores: 8, mem: 8, time: "02:00:00" }
bam_index:          { cores: 8, mem: 2, time: "02:00:00" }
bam_bw:             { cores: 8, mem: 2, time: "02:00:00" }
feature_counts:     { cores: 8, mem: 2, time: "01:00:00" }
edgeR:              { cores: 1, mem: 4, time: "04:00:00" }
edgeR_assemble:     { cores: 1, mem: 16, time: "04:00:00" }
juan:               { cores: 1, mem: 16, time: "01:00:00" }
```

Any rule without its own section falls back to `defaults`. Set `mapping.perform_mapping: False` if you're supplying your own BAM files and don't want splicekit to align FASTQs itself.

After setting up `splicekit.config` and `config.yaml`, run the pipeline as described in [Quick Start](quickstart.md).
