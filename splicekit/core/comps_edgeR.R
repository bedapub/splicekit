library("edgeR")
library("data.table")

args = commandArgs(trailingOnly = T); # trailingOnly: only take parameters after --args

input_folder = args[1]
feature_type = args[2]
control_name = args[3]
test_name = args[4]
comparison_name = args[5]
sample_membership = args[6]
sample_membership = unlist(strsplit(sample_membership, ","))
offset = 8
if (length(args) >= 7) {
    filter_low = args[7]
} else {
    filter_low = "filter_low" # filter lowly expressed features by default
}

right_offset = 0
#if (feature_type=="exons") { right_offset = 10 }
#if (feature_type=="junctions") { right_offset = 11 }

input_fname = paste(input_folder, "/data/samples_", feature_type, "_counts.tab.gz", sep="")
gx <- fread(input_fname, colClasses=list(character=1:1))
genes <- gx[, -c((offset+1):(ncol(gx)-right_offset)), with=FALSE]
gxcounts <- gx[, -c(1:offset), with=FALSE] # first 8 columns are feature infos
gxcounts <- gxcounts[, c(1:(ncol(gxcounts)-right_offset)), with=FALSE] # last 7 (or 10) columns are pji, delta_pji (PSI, delta_PSI)

group <- factor(sample_membership)
design <- model.matrix(~0 + group)
contrast <- makeContrasts(mytest = get(paste("group", test_name, sep="")) - get(paste("group", control_name, sep="")), levels = design)

dge = DGEList(counts=gxcounts, group=group, genes=genes)

# filter out features with low read counts
if (filter_low=="filter_low") { 
    keep = filterByExpr(dge, group=group)
    dge <- dge[keep, , keep.lib.sizes=FALSE]
}

dge = calcNormFactors(dge)
dge <- estimateDisp(dge, design, robust=TRUE)

fit <- glmQLFit(dge, design, robust=TRUE)

qlf <- glmQLFTest(fit, contrast=contrast)
d = topTags(qlf, n=Inf)
d = d$table
d["pi_value"] = -log10(d["PValue"]) * d["logFC"] # add pi value
output_fname = paste(input_folder, "/results/edgeR/", feature_type, "/", comparison_name, "_difffeature.tab.gz", sep="")
write.table(d, file=gzfile(output_fname), sep="\t", row.names=FALSE, quote=FALSE)

if (feature_type!="genes") {
    sp <- diffSpliceDGE(fit, geneid="gene_id", exonid="feature_id", contrast=contrast)
    d = topSpliceDGE(sp, test="exon", n=Inf)
    d["pi_value"] = -log10(d["P.Value"]) * d["logFC"] # add pi value
    output_fname = paste(input_folder, "/results/edgeR/", feature_type, "/", comparison_name, "_altsplice.tab.gz", sep="")
    write.table(d, file=gzfile(output_fname), sep="\t", row.names=FALSE, quote=FALSE)
}

# save dispersion plot
output_fname = paste(input_folder, "/results/edgeR/dispersion/", feature_type, "_dispersion.png", sep="")
png(output_fname)
plotBCV(dge)
dev.off()