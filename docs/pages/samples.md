# Sample annotation

The first step of the analysis, `splicekit annotation`, loads samples from `samples.tab` and builds the test-vs-control comparisons that every later step (features, edgeR, motifs, ...) works from.

!!! note
    `samples.tab` is **TAB delimited**. Lines starting with `#` are treated as comments.

```
sample_id	treatment_id
sample1	control
sample2	control
sample3	test1
sample4	test1
sample5	test2
sample6	test2
```

splicekit expects a BAM file per sample, resolved either from `bam_path` (as `{bam_path}/{sample_id}.bam`) or from a per-sample path in a `bam_column` column of `samples.tab` — see [Configuration](configuration.md#genome-parameters). The `sample_column`, `treatment_column`, `control_name`, `separate_column` and `group_column` parameters in `splicekit.config` control how `samples.tab` is read and how comparisons are formed.

## Comparisons

Each non-control treatment (which may have several replicate samples) is compared against the samples matching `control_name`. For the example above, this produces `annotation/comparisons.tab`:

```
comparison       compound_samples             DMSO_samples
test1_control    sample3_test1,sample4_test1  sample1_control,sample2_control
test2_control    sample5_test2,sample6_test2  sample1_control,sample2_control
```

In addition, `splicekit annotation` creates the processing shell scripts and (on `"cluster"`/`"SLURM"` platforms) cluster job files under `jobs/*`. An example cluster job file:

```bash
#!/bin/bash
#BSUB -J edgeR_junctions_sample1                  # job name
#BSUB -n 4                                        # number of tasks
#BSUB -R "span[hosts=1]"                          # 1 host
#BSUB -q short                                    # select queue
#BSUB -o logs_edgeR_junctions/sample1_control.out # output file
#BSUB -e logs_edgeR_junctions/sample1_control.err # error file

ml R
R --no-save --args splicekit comparison_junctions_data junctions control test ... < comps_edgeR.R
```

When running under Snakemake instead, job submission and resourcing come from `config.yaml` and `run_snakemake_slurm.sh` — see [Configuration](configuration.md#configyaml).
