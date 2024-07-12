library(DEXSeq)
library("data.table")
library(BiocParallel)
register(MulticoreParam(6))

args = commandArgs(trailingOnly = T); # trailingOnly: only take parameters after --args

input_folder = args[1]
feature_type = args[2]
control_name = args[3]
test_name = args[4]
comparison_name = args[5]
sample_membership = args[6]
sample_membership = unlist(strsplit(sample_membership, ","))
samples = args[7]
samples = unlist(strsplit(samples, ","))

sample_counts <- as.list(paste("results/dexseq/counts/", samples, "_dexseq_counts.txt.clean", sep = ""))

countFiles <- sample_counts
sampleTable <- data.frame(
    row.names = samples,
    condition = sample_membership,
    libType = as.list(rep("paired-end", length(sample_membership)))
)

dxd <- DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData = sampleTable,
    design = ~ sample + exon + condition:exon,
    flattenedfile = "dexseq.gff"
)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM = MulticoreParam(6))
dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition")
result <- DEXSeqResults(dxd)

result_df <- as.data.frame(result)
write.table(result_df, file=gzfile("results/dexseq/dexseq_results.tab.gz"), sep="\t", quote=FALSE, row.names=FALSE)