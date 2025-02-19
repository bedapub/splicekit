#!/bin/bash
snakemake --executor cluster-generic --cluster-generic-submit-cmd "sbatch --job-name=splicekit-{rule}-{wildcards} --output={log}.out --error={log}.err --parsable -n {resources.cores} --mem-per-cpu={resources.mem}G --time={resources.time}" --cluster-generic-cancel-cmd "scancel" --jobs 100 "$@"
