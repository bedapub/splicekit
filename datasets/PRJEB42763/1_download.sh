# downloads fastq files for this dataset
# requires fastq-dump

# Usually, it can be loaded on a cluster node with:
# ml .testing; ml SRA-Toolkit/2.11.0-gompi-2020a

mkdir fastq

runs="ERR5296631
ERR5296641
ERR5296646
ERR5296647
ERR5296827
ERR5296849
ERR5296853
ERR5296899
ERR5300363
ERR5301820
ERR5302785
ERR5303450
ERR5303456
ERR5303461
ERR5303484
ERR5303525
ERR5303527
ERR5303529
ERR5305084"

cd fastq

for rid in $runs
do
    echo "downloading $rid"
    fastq-dump --split-3 -A $rid
    pigz -f ${rid}_1.fastq
    pigz -f ${rid}_2.fastq
done
