# maps fastq files to the homo_sapiens reference genome
# requires STAR

# download and prepare the homo_sapiens latest version genome from Ensembl
pybio genome homo_sapiens

mkdir bam

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

for rid in $runs
do
    echo "mapping $rid"
    pybio star homo_sapiens fastq/${rid}_1.fastq.gz fastq/${rid}_2.fastq.gz bam/${rid}
done
