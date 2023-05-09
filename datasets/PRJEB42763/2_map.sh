# maps fastq files to the homo_sapiens reference genome
# requires STAR

# download and prepare the homo_sapiens latest version genome from Ensembl
pybio genome homo_sapiens

mkdir bam

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

for rid in $runs
do
    echo "mapping $rid"
    pybio star homo_sapiens fastq/${rid}_1.fastq.gz fastq/${rid}_2.fastq.gz bam/${rid}
done
