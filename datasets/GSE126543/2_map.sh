# maps fastq files to the homo_sapiens reference genome
# requires STAR

# download and prepare the homo_sapiens latest version genome from Ensembl
pybio genome homo_sapiens

mkdir bam

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

for rid in $runs
do
    echo "mapping $rid"
    pybio star homo_sapiens fastq/${rid}_1.fastq.gz fastq/${rid}_2.fastq.gz bam/${rid}
done
