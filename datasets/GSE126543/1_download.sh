# downloads fastq files for this dataset
# requires fastq-dump

# Usually, it can be loaded on a cluster node with:
# ml .testing; ml SRA-Toolkit/2.11.0-gompi-2020a

mkdir fastq

runs="SRR8571937
SRR8571938
SRR8571939
SRR8571940
SRR8571941
SRR8571942
SRR8571943
SRR8571944
SRR8571945
SRR8571946
SRR8571947
SRR8571948
SRR8571949
SRR8571950
SRR8571951
SRR8571952"

cd fastq

for rid in $runs
do
    echo "downloading $rid"
    fastq-dump --split-3 -A $rid
    gzip -f ${rid}_1.fastq
    gzip -f ${rid}_2.fastq
done
