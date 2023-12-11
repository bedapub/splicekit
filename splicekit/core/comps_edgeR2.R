library("edgeR")
library("data.table")

args = commandArgs(trailingOnly = T); # trailingOnly: only take parameters after --args

input_folder = args[1]
atype = args[2]
control_name = args[3]
test_name = args[4]
comp_name = args[5]
samples_control = args[6]
samples_control = unlist(strsplit(samples_control, ","))
samples_test = args[7]
samples_test = unlist(strsplit(samples_test, ","))
num_control = as.numeric(length(samples_control))
num_test = as.numeric(length(samples_test))
offset = 8
if (length(args) >= 9) {
    filter_low = args[9]
} else {
    filter_low = "filter_low" # filter our lowly expressed features by default
}

samples_all = c(samples_control, samples_test)
samples_all_plain <- gsub("_kd|_control", "", samples_all)

right_offset = 0
#if (atype=="exons") { right_offset = 10 }
#if (atype=="junctions") { right_offset = 11 }

input_fname = paste(input_folder, "/data/samples_", atype, "_counts.tab.gz", sep="")
gx <- fread(input_fname, colClasses=list(character=1:1))
genes <- gx[, -c((offset+1):(ncol(gx)-right_offset)), with=FALSE]
gxcounts <- gx[, -c(1:offset), with=FALSE] # first 8 columns are feature infos
gxcounts <- gxcounts[, c(1:(ncol(gxcounts)-right_offset)), with=FALSE] # last 7 (or 10) columns are pji, delta_pji (PSI, delta_PSI)

y = DGEList(counts=gxcounts, genes=genes)
y = calcNormFactors(y)
#y <- estimateDisp(y, robust=TRUE)
#print(y$common.dispersion)

y_subset = y[, colnames(y) %in% c(samples_control, samples_test)]

listA = colnames(y_subset)
listB = samples_control
listC = samples_test
list_member = character(length(listA))
for (i in seq_along(listA)) {
  if (listA[i] %in% listB) {
    list_member[i] = "control"
  } else if (listA[i] %in% listC) {
    list_member[i] = "test"
  } else {
    list_member[i] = "None"
  }
}

# with intercept
#group <- factor(list_member, levels = c("control", "test"), labels = c("control", "test"))
#design <- model.matrix(~group)
#colnames(design) <- c("intercept", "test")
#contrast <- makeContrasts(mytest = test, levels = design)
#contrast

# without intercept
group <- factor(list_member, levels = c("control", "test"), labels = c("control", "test"))
design <- model.matrix(~0 + group)
colnames(design) <- c("control", "test")
contrast <- makeContrasts(mytest = test - control, levels = design)
contrast

# filter out features with low read counts
if (filter_low=="filter_low") { 
    keep <- filterByExpr(y_subset, group=group, min.count = 5, min.total.count = 10, large.n = 3) 
    table(keep)
    y_subset <- y_subset[keep,]
}

y_subset <- estimateDisp(y_subset, design, robust=TRUE)

fit <- glmQLFit(y_subset, design, robust=TRUE)

qlf <- glmQLFTest(fit, contrast=contrast)
d = topTags(qlf, n=Inf)
d = d$table
d["pi_value"] = -log10(d["PValue"]) * d["logFC"] # add pi value
output_fname = paste(input_folder, "/results/results_edgeR2_", atype, "/", comp_name, "_difffeature.tab.gz", sep="")
write.table(d, file=gzfile(output_fname), sep="\t", row.names=FALSE, quote=FALSE)

sp <- diffSpliceDGE(fit, geneid="gene_id", exonid="feature_id", contrast=contrast)
d = topSpliceDGE(sp, test="exon", n=Inf)
d["pi_value"] = -log10(d["P.Value"]) * d["logFC"] # add pi value
output_fname = paste(input_folder, "/results/results_edgeR2_", atype, "/", comp_name, "_altsplice.tab.gz", sep="")
write.table(d, file=gzfile(output_fname), sep="\t", row.names=FALSE, quote=FALSE)
