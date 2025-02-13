mkdir fastq

runs="SRR22909639
SRR22909638
SRR22909637
SRR22909636
SRR22909635
SRR22909634
SRR22909633
SRR22909632
SRR22909631
SRR22909630
SRR22909629
SRR22909628
SRR22909627
SRR22909626
SRR22909625"

cd fastq

for rid in $runs
do
    echo "downloading $rid"
    fastq-dump --split-3 -A $rid
    gzip -f ${rid}_1.fastq
    gzip -f ${rid}_2.fastq
done
