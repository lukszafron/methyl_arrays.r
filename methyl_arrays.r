#! /bin/env Rscript

options(max.print = 10000)
cat("Program path:", unlist(strsplit(grep(commandArgs(), pattern = "file=", value = T), split = "="))[2], "\n")

arguments <- commandArgs(trailingOnly = T)
# arguments <- c('/shares/Microarrays/OvCa', '/shares/Microarrays/OvCa/OvCa.microarrays.final.csv', '/workspace/lukasz/Microarrays', "TRUE", "60", "Sample_Name", "Sample_Group+Sample_Source", "/workspace/lukasz/Microarrays/OvCa.microarrays.final/microarrays.final.sig_CpGs.txt", "FALSE")
# arguments <- c("/shares/Microarrays/Macice/", "/shares/Microarrays/Macice/Macice.csv", "/workspace/lukasz/Microarrays", "TRUE", "60", "Sample_Name", "Sample_Group", "NA", "FALSE")
if(length(arguments) != 9) {cat(paste("The number of provided arguments is incorrect:", 
                                       length(arguments), "instead of 9.
                                       The arguments should be placed in the following order:
                                       1 - a directory where microarray data are stored,
                                       2 - a csv file with microarray data description,
                                       3 - a directory where the analysis results should be saved,
                                       4 - a boolean value (TRUE/FALSE) determining if the data binarization step should be performed,
                                       5 - a number of CPU threads to be used,
                                       6 - a variable in the aforementioned csv file containing sample names,
                                       7 - independent factor variables in the aforementioned csv file to be used, separated with '+', starting with a grouping variable,
                                       8 - a path to the optional txt file containing the list of CpG sites (one per a line) for which the methylation status visualization is to be performed (NA if none),
                                       9 - a boolean value (TRUE/FALSE) indicating if the cohort of patients contains individuals of both sexes.\n"))
  stop("Incorrect number of arguments.")}
arguments.backup <- arguments

read.args <- function() {
dataDirectory <<- arguments[1]
csvfile <<- arguments[2]
workspace <<- arguments[3]
illuminaDataDir <<- file.path(workspace, "Illumina_data")
dataBinarization <<- as.logical(arguments[4])
threads <<- as.numeric(arguments[5])
snames <<- arguments[6]
ind.factors <<- unlist(strsplit(arguments[7], split = "\\+"))
CpGs.signature.file <<- arguments[8]

if(CpGs.signature.file != "NA") {
sel.CpGs <<- readLines(con = file(CpGs.signature.file))
close(file(CpGs.signature.file))} else {
sel.CpGs <<- "NA"
}
both.sexes <<- as.logical(arguments[9])
}
read.args()

genome.ver <- "hg19" # Do not modify this value
methyl.type <- "EPIC"
suffix <- paste("dataDir", basename(dataDirectory), "clinData", basename(csvfile), "dataBin", dataBinarization, "methyl.type", methyl.type, "genome", genome.ver, "both.sexes", both.sexes, sep = ".")

workdir <- sub(basename(csvfile), pattern = "\\.[^\\.]*$", replacement = "")
workspace <- file.path(workspace, workdir)

cross_hybrid.probes.txt <- file.path(illuminaDataDir, "Illumina_EPIC_potentially_cross-hybridising_CpG_targeting_probes_PMID:27330998.txt")
refseq.bed <- file.path(illuminaDataDir, "hg19_NCBI_RefSeq.bed")
dhs.bed <- file.path(illuminaDataDir, "hg19_wgEncodeRegDnaseClusteredV3.bed")
tfbs.bed <- file.path(illuminaDataDir, "/hg19_wgEncodeRegTfbsClusteredV3.bed")
chain <- file.path(illuminaDataDir, "hg19ToHg38.over.chain")

if(dataBinarization) {
  log2ratios.list <- c("log2ratios", "log2ratios.bin")} else {
  log2ratios.list <- c("log2ratios")
}

dir.create(workspace, recursive = T)
setwd(workspace)

library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(missMethyl)
library(Gviz)
library(DMRcate)
library(stringr)
library(tidyverse)
library(ggfortify)
library(pheatmap)
library(ggplot2)
library(ggpmisc)
library(pals)
library(factoextra)
library(foreach)
library(doMC)
registerDoMC(threads)
library(openxlsx)
library(pdftools)
library(UpSetR)
library(rtracklayer)
library(AnnotationHub)
library(GenomicFeatures)
library(data.table)
library(annotatr)
library(plyr)
library(dplyr)
library(tibble)
library(stringi)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)

if(!file.exists(paste("Methylation analysis results", suffix, "RData", sep = "."))) {
  
annots <- c("hg19_basicgenes", "hg19_genes_intergenic", "hg19_genes_intronexonboundaries", "hg19_lncrna_gencode", "hg19_genes_firstexons", "hg19_genes_cds", "hg19_enhancers_fantom")
annotations <- build_annotations(genome = genome.ver, annotations = annots)

dataTable <- fread(csvfile)
sink(paste("Analysis summary", suffix, "txt", sep = "."))
cat("Total number of samples:", nrow(dataTable), "\n\n")
print(with(dataTable, by(dataTable, INDICES = lapply(dataTable[,..ind.factors], as.vector), FUN = function(x) {l1 <-list(nrow(x), x[[snames]]); names(l1) <- c("Number of samples", "List of samples"); l1})))
sink()

targets <- read.metharray.sheet(dataDirectory, pattern = basename(csvfile)) # Define targets
ah <- AnnotationHub(ask = F)
setAnnotationHubOption("max_downloads", threads)
epiFiles <- query(ah, "EpigenomeRoadMap")
GencodeV10 <- epiFiles[["AH49010"]] # Get all gene and transcript annotations for the hg19 genome from Gencode
ch19_38 <- import.chain(con = chain)

EPICanno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(EPICanno)

for(i in ind.factors) {targets[[i]] <- factor(targets[[i]])}

cat("Reading microarray data into the R rgSet object...\n")
rgSet <- read.metharray.exp(targets=targets, force = T) # Read data into R
targets$ID <- paste(targets[[ind.factors[1]]], targets[[snames]], sep = ".")
sampleNames(rgSet) <- targets$ID
print(rgSet)

# calculate the detection p-values
cat("Calculating the detection p-values...\n")
detP <- detectionP(rgSet)

# examine mean detection p-values across all samples to identify any failed samples
colMeans.detP <- data.frame(colMeans(detP))
colMeans.detP <- setNames(colMeans.detP, nm = "colMeans")
colMeans.detP$SD <- colSds(detP)
colMeans.detP$SE <- colMeans.detP$SD/sqrt(nrow(detP))
t = qt((1-0.05)/2 + .5, nrow(detP)-1) # Value of the Student's t distribution for alpha = 0.05, tends to be 1.96 if the sample size is big enough.
colMeans.detP$CI95 <- t*colMeans.detP$SE

qcReport(rgSet, sampNames=targets$ID, sampGroups=targets[[ind.factors[1]]], pdf=paste("qcReport", suffix, "pdf", sep = "."))

if(nrow(colMeans.detP) < 56) {df.width <- 7} else {df.width <- nrow(colMeans.detP)/8}
pdf(paste("Microarray_quality_control", suffix,"pdf", sep = "."), width = df.width)
quality.plot <- ggplot(colMeans.detP,aes(x=reorder(rownames(colMeans.detP), as.numeric(targets[[ind.factors[1]]])), y=colMeans, color=targets[[ind.factors[1]]])) + 
  geom_point(stat="identity") + scale_y_log10(breaks = c(0.001,0.01, 0.05, 0.1)) + 
  scale_color_brewer(palette = "Set1") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_errorbar(aes(y = colMeans, ymin=colMeans-CI95, ymax=colMeans+CI95)) + 
  geom_hline(yintercept=0.05, linetype="dashed", color = "red") + 
  xlab("Samples") + ylab("Means of detection p-values") + 
  labs(title = "Microarray quality assessment", color = ind.factors[1], subtitle = "(with 95% confidence intervals)") + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
print(quality.plot)

detP.tibble <- as_tibble(detP)
detP.centiles <- foreach(i = colnames(detP.tibble), .combine = cbind) %dopar% {temp <- quantile(as.data.frame(detP.tibble %>% dplyr::select(i))[,1], seq(0,1,0.05))}

colnames(detP.centiles) <- make.names(sub(colnames(detP.tibble),pattern = "_df", replacement = ""))
rownames(detP.centiles) <- as.numeric(sub(rownames(detP.centiles), pattern = "%", replacement = ""))
detP.centiles.df <- as.data.frame(detP.centiles)

keep <- detP.centiles.df[rownames(detP.centiles.df) == "85",] < 0.05
if(!all(foreach(i = rgSet[[snames]], j = colnames(keep), .combine = c) %do% {grepl(j, pattern = i)})) {stop("The order of samples in the rgSet and keep objects do not match.")}

keep <- as.vector(keep)
rgSet <- rgSet[, keep] # The poor-quality samples are being filtered out.
targets <- targets[keep, ] # Remove filtered out samples.
detP <- detP[, keep] # Remove poor quality samples from detection p-value table
if (!all(sampleNames(rgSet) == targets[["ID"]])) {stop("The sample names do not mach after filtering.")}
if (!all(sampleNames(rgSet) == colnames(detP))) {stop("The sample names do not mach after filtering.")}

sink(paste("Analysis summary", suffix, "txt", sep = "."), append = T)
cat("\nSamples filtered out due to their poor quality:\n\n")
print(dataTable[!dataTable[[snames]] %in% targets[[snames]],c(1,3,4)])
sink()

if(ncol(detP.centiles.df) > 16) {
df.nos <- ceiling(ncol(detP.centiles.df)/16)} else {df.nos <- 1}
for(i in seq(1, df.nos)) { if(i*16 <= ncol(detP.centiles.df)) {tmp.df <- detP.centiles.df[,seq(i*16-15,i*16)]} else{
  tmp.df <- detP.centiles.df[,seq(i*16-15, ncol(detP.centiles.df)), drop = F]}
  assign(paste("qdf", i, sep = "."), value = tmp.df)}
for (j in ls(pattern = "^qdf.[0-9]+")) {
dflist <- colnames(get(j))
quality.plot2.tmp <- ggplot(get(j)) + 
  sapply(dflist, function(i) {geom_line(aes_string(x = as.numeric(rownames(get(j))), y = i, color = 'i', group = 'i'))}) + 
  scale_color_manual(values=as.vector(cols25())) +
  scale_y_log10() + scale_y_continuous(breaks = c(0.01, 0.05, 0.1, 0.5, 1)) +
  scale_x_continuous(breaks = seq(0,100,5)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red") +
  xlab("Percentile of methylation sites") + ylab("Detection p-values") +
  labs(title = "Microarray quality assessment 2", color = "Samples") +
  theme(plot.title = element_text(hjust = 0.5))
  print(quality.plot2.tmp)}
dev.off()

rm(list = grep(ls(), pattern = "^qdf.[0-9]+", value = TRUE))

cat("Performing the data normalization step...\n")
mSetSq <- preprocessFunnorm(rgSet = rgSet) # Data normalization step. As a results, we get the GenomicRatioSet object, with M-values, being the log2 values of Meth/Unmeth fluorescence after normalization. To retrieve this ratios, we use the getM() function.
# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet) # We get the MethylSet object. To retrive the metyl values, use the getMeth() function, and to get the unmetyl. values, use the getUnmeth() function.

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
if(!all(featureNames(mSetSq) == rownames(detP))) {stop("The feature names do not match.")}

sink(paste("Analysis summary", suffix, "txt", sep = "."), append = T)
cat("\n\nThe total number of probes available in the microarray was ", nrow(mSetSq),".\n", sep = "")
sink()
# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.05) == ncol(mSetSq) 
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

filter.summary <- function(x, filter) {cat("The number of probes passing the ", filter, " filter was: ", x, " (", round(x/nrow(mSetSq)*100,2),"%).\n",sep = "")}

sink(paste("Analysis summary", suffix, "txt", sep = "."), append = T)
filter.summary(nrow(mSetSqFlt), "'detP <0.05 for all the samples'")
sink()

if(both.sexes) {
# # if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(mSetSqFlt) %in% EPICanno$Name[EPICanno$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
}
# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

sink(paste("Analysis summary", suffix, "txt", sep = "."), append = T)
filter.summary(nrow(mSetSqFlt), "'SNPs at CpG site'")
sink()

# exclude cross reactive probes
cross_hybrid.probes <- readLines(con = file(cross_hybrid.probes.txt))
on.exit(close(con = file(cross_hybrid.probes.txt)))
keep <- !featureNames(mSetSqFlt) %in% cross_hybrid.probes
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt

sink(paste("Analysis summary", suffix, "txt", sep = "."), append = T)
filter.summary(nrow(mSetSqFlt), "'cross reactive probes'")
sink()

log2ratios <- getM(mSetSqFlt) # Get log2 values of Methyl/Unmethyl ratios.
Beta.values <- getBeta(mSetSqFlt) # Get Beta values (i.e. Methyl/(Methyl+Unmethyl))
log2ratios[is.infinite(log2ratios)] <- NA
#Delete all probes with at least one NA value.
Beta.values <- Beta.values[rowSums(is.na(log2ratios)) == 0,]
log2ratios <- log2ratios[rowSums(is.na(log2ratios)) == 0,]

sink(paste("Analysis summary", suffix, "txt", sep = "."), append = T)
filter.summary(nrow(log2ratios), "'infinite M value'")
sink()

# Binarize the log2ratios
log2ratios.bin <- log2ratios
log2ratios.bin[log2ratios.bin >= 0] <- 1
log2ratios.bin[log2ratios.bin < 0] <- -1

# Binarize the Beta values
Beta.values.bin <- Beta.values
Beta.values.bin[Beta.values.bin >= 0.5] <- 1
Beta.values.bin[Beta.values.bin < 0.5] <- 0

#Draw density plots of Beta-values
cat("Drawing density plots of beta values...\n")
# visualise what the data looks like before and after normalisation and filtering
pdf(paste("Density_plots", suffix, "pdf", sep = "."), height = 10, width = 7)
par(mfrow=c(3,1))
densityPlot(getBeta(mSetRaw), sampGroups=targets[[ind.factors[1]]],
            main="Raw dataset", xlab = "Beta-values", legend=FALSE, pal = cols25())
legend("topright", legend = levels(factor(targets[[ind.factors[1]]])), 
       text.col=cols25(), cex = 1, bg = NULL)
abline(v = 0.5, col = "red", lty = "dashed")

densityPlot(getBeta(mSetSqFlt), sampGroups=targets[[ind.factors[1]]],
            main="Normalized and filtered dataset", xlab = "Beta-values", legend=FALSE, pal = cols25())
legend("topright", legend = levels(factor(targets[[ind.factors[1]]])), 
       text.col=cols25(), cex = 1, bg = NULL)
abline(v = 0.5, col = "red", lty = "dashed")

densityPlot(Beta.values.bin, sampGroups=targets[[ind.factors[1]]],
            main="Normalized and filtered dataset after binarization",
            xlab = "Beta-values", legend=FALSE, pal = cols25())
legend("topright", legend = levels(factor(targets[[ind.factors[1]]])), 
      text.col=cols25(), cex = 1, bg = NULL)
abline(v = 0.5, col = "red", lty = "dashed")

dev.off()

#Calculate distances between samples
cat("Plotting heatmaps portraying distances between samples...\n")
if(dim(log2ratios)[2]>35) {width <- ceiling(dim(log2ratios)[2]/5)} else {width <- 7}
height <- width
pdf(paste("Heatmap_sample_distances", suffix, "pdf", sep = "."), width = width, height = height)

for(l2ratios in log2ratios.list){
sampleDists <- dist(t(get(l2ratios)))
sampleDistMatrix <- as.matrix(sampleDists)

# Plot heatmaps portraying distances between samples

anno <- targets %>% dplyr::select(c("ID", ind.factors))
anno <- as.data.frame(lapply(anno, as.factor))
rownames(anno) <- anno$ID
anno[["ID"]] <- NULL
levels.no <- function(x) {length(levels(x))}
cols <- brewer.paired(12)
col.pals <- foreach(col.pal = sapply(anno, levels.no)) %do% {sel.cols <- cols[1:col.pal]; cols <- rev(cols); sel.cols}
names(col.pals) <- colnames(anno)
invisible(foreach(var.no = seq(1,length(col.pals)), .combine = c) %do% {names(col.pals[[var.no]]) <- as.vector(sapply(anno[var.no], levels))})

if(!all(rownames(sampleDistMatrix) == rownames(anno))) {stop("The sample names do not match.")}

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = anno,
         annotation_colors = col.pals,
         cellwidth = 10, cellheight = 10, main = 
         if(l2ratios == "log2ratios"){
            main.title <- "Original sample distances"; main.title} else {
            main.title <- "Binarized sample distances"; main.title
          }
          )
}
dev.off()

## PCA analysis
cat("Performing principal component analysis...\n")
pdf(paste("PCA_analysis_plots", suffix, "pdf", sep = "."))
for(l2ratios in log2ratios.list){
log2ratios.t <- t(get(l2ratios))
log2ratios.t.pca.data <- prcomp(log2ratios.t)

# Show the eigenvalues
print(fviz_eig(log2ratios.t.pca.data, addlabels = T, main = paste("Scree plot for", if(l2ratios == "log2ratios"){"M-values before binarization"} else {"M-values after binarization"})))

# Show the PCA plot for both groups and sources
assign(paste(l2ratios, "pca.data", sep = "."), value = log2ratios.t.pca.data, envir = .GlobalEnv)
log2ratios.t.df <- as.data.frame(log2ratios.t)
log2ratios.t.df.merged <- merge(x = log2ratios.t.df, y = anno, by.x = 0, by.y = 0, sort = F)
rownames(log2ratios.t.df.merged) <- log2ratios.t.df.merged[["Row.names"]]
log2ratios.t.df.merged <- log2ratios.t.df.merged[,-1]
log2ratios.t.df.merged <- log2ratios.t.df.merged[match(rownames(log2ratios.t.pca.data$x), rownames(log2ratios.t.df.merged)), , drop = F]

ind.length <- sapply(anno, levels.no)
if(ncol(log2ratios.t.pca.data$x) > 1) {
  if(length(ind.length) >1) {
  if(ind.length[1] <= 6) {
    colors.n <- ind.length[2]
    ap1 <- autoplot(log2ratios.t.pca.data, data = log2ratios.t.df.merged, size = 3, shape = names(ind.length[1]), colour = names(ind.length[2]))} else
  if(ind.length[2] <= 6) {
    colors.n <- ind.length[1]
    ap1 <- autoplot(log2ratios.t.pca.data, data = log2ratios.t.df.merged, size = 3, shape = names(ind.length[2]), colour = names(ind.length[1]))} else {
      colors.n <- ind.length[1] * ind.length[2]
      ap1 <- autoplot(log2ratios.t.pca.data, data = log2ratios.t.df.merged, size = 3, colour = interaction(names(ind.length[1]), names(ind.length[2])))}
  } else {
    colors.n <- ind.length[1]
    ap1 <- autoplot(log2ratios.t.pca.data, data = log2ratios.t.df.merged, size = 3, colour = names(ind.length[1]))
  }
    if(l2ratios == "log2ratios"){
    pca.title = paste("PCA plot - original data set")} else {
      pca.title = paste("PCA plot - binarized data set")
    }
  ap1 <- ap1 + scale_color_manual(values = if(colors.n > 25) {rainbow(colors.n)} else {as.vector(cols25())}) + coord_fixed() + 
    labs(title = pca.title) +
    theme(plot.title = element_text(hjust = 0.5))
  try(print(ap1))
  } else {
    plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main="ERROR: Too few dimensions to draw a PCA plot.")
  }

# Show the PCA plot for specific samples
print(fviz_pca_ind(log2ratios.t.pca.data,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = F) + labs(title = paste("PCA with qualities of representation", if(l2ratios == "log2ratios"){"(M-values before binarization)"} else {"(M-values after binarization)"})))

# Results for Variables
res.var <- get_pca_var(log2ratios.t.pca.data)
# Get 10 the most contributing vars for dim1 and 10 for dim2.
sig.cgs <- c(names(head(sort(res.var$contrib[,1], decreasing = T), 10)), names(head(sort(res.var$contrib[,2], decreasing = T), 10)))
l2ratios.pca.sel <- prcomp(t(get(l2ratios)[sig.cgs,]))
anno.prcomp <- anno[match(rownames(l2ratios.pca.sel$x), rownames(anno)), , drop = F]

# Show the biplot for these 20 variables.
print(fviz_pca_biplot(l2ratios.pca.sel, repel = F,
                col.var = "#2E9FDF", # Variables color
                col.ind = anno.prcomp[[ind.factors[1]]]
) + labs(title = paste("PCA - Biplot for 20 the most contributing variables", if(l2ratios == "log2ratios"){"(M-values before binarization)"} else {"(M-values after binarization)"})))
}
dev.off()

# Subsetting the EPICanno annotation object
EPICannoSub <- EPICanno[match(rownames(log2ratios), rownames(EPICanno)),]

# Create genomic ranges object from the Illumina EPIC db file.
cpgData <- GRanges(seqnames=Rle(EPICanno$chr),
                   ranges=IRanges(start=EPICanno$pos, end=EPICanno$pos),
                   strand=Rle(EPICanno$strand))
names(cpgData@ranges) <- rownames(EPICanno)

# Add more detailed annotations with the use of the "annotatr" library

if(!file.exists(file.path(illuminaDataDir, "Illumina_Infinium_methyl_EPIC_array_hg19_extended_attributes.csv"))) {
  cat("Assigning additional attributes to the Illumina methylation data...\n")
  
  if(threads >30) {threads <- 30}
 registerDoMC(threads) 
  cpgData.annotated <- as.data.frame(rbindlist(foreach(cpg.no = seq(1, length(cpgData))) %dopar% {
    cpgData.sel <- cpgData[cpg.no]
    strand(cpgData.sel) <- "*"
    anno.cpg.no <- subsetByOverlaps(annotations, cpgData.sel)
    anno.gene <- setNames(as.data.frame(t(sapply(strsplit(gsub(gsub(gsub(gsub(paste(paste(paste0(unique(paste(as.data.frame(anno.cpg.no)[["type"]], 
                                                                                                              as.data.frame(anno.cpg.no)[["strand"]], sep = "(")), ")"), 
                                                                                          collapse = ";"), paste(unique(paste0(paste(as.data.frame(anno.cpg.no)[["symbol"]], 
                                                                                                                                     as.data.frame(anno.cpg.no)[["strand"]], 
                                                                                                                                     sep = "("), ")")), 
                                                                                                                 collapse = ";"), sep = "#"), 
                                                                              pattern = "\\(-\\)", replacement = "(m)"), 
                                                                         pattern = "\\(\\+\\)", replacement = "(p)"), 
                                                                    pattern = "\\(\\*\\)", replacement = "(b)"), 
                                                               pattern = "(hg19_|genes_|_gencode|_fantom)", replacement = ""), split = "#"), cbind))), nm = c("Annotations", "Genes"))
    anno.gene <- as.data.frame(sapply(anno.gene, FUN = function(x) {sub(x = x, pattern = "^\\)$", replacement = "NA")}, simplify = F))
    
    if(length(anno.cpg.no) == 0) {anno.gene <- anno.gene %>% dplyr::mutate(across(.fns = function(x) {paste0(x, "(b)")}))}
    
    return(as.data.frame(cpgData[cpg.no]) %>% add_column(anno.gene, .after = 5))
  }, use.names = T))
  
  threads <- as.numeric(arguments[5])
  registerDoMC(threads)
  
  cpgData.annotated <- cpgData.annotated[order(cpgData.annotated$seqnames, cpgData.annotated$start),]
  unique.genes <- sort(unique(unlist(strsplit(cpgData.annotated$Genes, split = ";"))))
  ug <- unique(sub(unique.genes, pattern = "\\(.*\\)$", replacement = ""))
  ug <- grep(ug, pattern = "^NA$", invert = T, value = T)
  ambiguous.genes.no <- foreach(g = ug, .combine = c) %dopar% {v1 <- grep(unique.genes, pattern = paste0("(^|;)",g,"\\(.*\\)")); if(length(v1) > 1) {v1}}
  if(!is.null(ambiguous.genes.no)) {
  ambiguous.genes.pattern <- gsub(paste0("(", paste(sub(unique.genes[ambiguous.genes.no], pattern = "(.*)", replacement = "(^|;)\\1(;|$)"), collapse = "|"), ")"), pattern = "\\(([bmp])\\)", replacement = "\\\\(\\1\\\\\\)")
  cpgData.annotated$Genes <- gsub(gsub(cpgData.annotated$Genes, pattern = ambiguous.genes.pattern, replacement = ";NA(b);"), pattern = ambiguous.genes.pattern, replacement = ";NA(b);")
  }
  cpgData.annotated$Genes <- gsub(cpgData.annotated$Genes, pattern = "((^;)|(;$))", replacement = "")
  
  cpgData.annotated$Genes.stranded <- foreach(index = seq(nrow(cpgData.annotated)), .combine = c) %dopar% {
    if(cpgData.annotated[["strand"]][index] == "+") {strand.abbr <- "\\(p\\)"} else 
      if(cpgData.annotated[["strand"]][index] == "-") {strand.abbr <- "\\(m\\)"} else 
      {stop("The strand identification has failed.")}
    unlist(sapply(strsplit(cpgData.annotated[["Genes"]][index], split = ";"), 
                  FUN = function(x) {if(any(grepl(x, pattern = strand.abbr))) {vec1 <- gsub(grep(x, pattern = strand.abbr, value = T), pattern = strand.abbr, replacement = "")
                  vec1 <- vec1[vec1 != "NA"]
                  if(length(vec1) == 0) {vec1 <- NA} else
                    if(length(vec1) > 1){
                      vec1 <- paste(vec1, collapse = "&")}
                  vec1} else {NA}}))
  }
  
  cat("Converting gene symbols' aliases to their primary names for all CpGs...\n")
  
  gene.list.conv <- NULL
  for(i in cpgData.annotated$Genes.stranded) {
    gene.list.conv <- append(gene.list.conv, values = paste(unique(foreach(gene = unlist(strsplit(i, split = "&")), .combine = c) %do% {
      new.symbol <- alias2Symbol(gene, species = "Hs")
      if(length(new.symbol) == 0) {
        gene} else {new.symbol}}), collapse = "&"))
  }
  
  cpgData.annotated$Genes.stranded <- gene.list.conv
  
  cat("Mapping gene symbols to ENSEMBL ids for all CpGs...\n")
  
  null.output <- file("/dev/null")
  sink(null.output, type = "output")
  sink(null.output, type = "message")
  
  cpgData.annotated$Genes.stranded.ENSEMBL <- foreach(genes = strsplit(cpgData.annotated$Genes.stranded, split = "&"), .combine = c) %do% {
    if(any(!is.na(genes))) {try.res <- try(
      expr = {mapped.genes <- mapIds(x = org.Hs.eg.db, keys = genes, column = "ENSEMBL", 
                                     keytype = "SYMBOL", multiVals = "first")
      paste(mapped.genes[!is.na(mapped.genes)], collapse = "&")})
    if(class(try.res) == "try-error") {
      if(grepl(try.res[1], pattern = "None of the keys entered are valid keys for 'SYMBOL'")) {
        NA
      } else {stop("The mapIds function failed for genes: ", paste(genes, collapse = ","), "\n", try.res[1])}
    } else if(class(try.res) == "character") {
      if(nchar(try.res) >0) {try.res} else {NA}
    } else {stop("An unidentified error in the gene mapping process occurred.")}
    } else {NA}
  }
  sink()
  sink()
  
  cat("Mapping gene symbols to ENTREZ ids for all CpGs...\n", sep = "")
  
  sink(null.output, type = "output")
  sink(null.output, type = "message")
  
  cpgData.annotated$Genes.stranded.ENTREZID <- foreach(genes = strsplit(cpgData.annotated$Genes.stranded, split = "&"), .combine = c) %do% {
    if(any(!is.na(genes))) {try.res <- try(
      expr = {mapped.genes <- mapIds(x = org.Hs.eg.db, keys = genes, column = "ENTREZID", 
                                     keytype = "SYMBOL", multiVals = "first")
      paste(mapped.genes[!is.na(mapped.genes)], collapse = "&")})
    if(class(try.res) == "try-error") {
      if(grepl(try.res[1], pattern = "None of the keys entered are valid keys for 'SYMBOL'")) {
        NA
      } else {stop("The mapIds function failed for genes: ", paste(genes, collapse = ","), "\n", try.res[1])}
    } else if(class(try.res) == "character") {
      if(nchar(try.res) >0) {try.res} else {NA}
    } else {stop("An unidentified error in the gene mapping process occurred.")}
    } else {NA}
  }
  sink()
  sink()
  
  write.table(x = cpgData.annotated, row.names = T, col.names = NA, file = file.path(illuminaDataDir, "Illumina_Infinium_methyl_EPIC_array_hg19_extended_attributes.csv"), sep = ";")
} else {
  cat("Loading the pre-existing table with extended gene attributes...\n", sep = "")
  cpgData.annotated <- read.table(file.path(illuminaDataDir, "Illumina_Infinium_methyl_EPIC_array_hg19_extended_attributes.csv"), sep = ";", header = T, row.names = 1)
}

unique.genes <- sort(unique(unlist(strsplit(cpgData.annotated$Genes, split = ";"))))
unique.genes <- grep(unique.genes, pattern = "(^|;)NA\\(", invert = T, value = T)

# Create a design matrix
design <- with(targets,model.matrix(eval(parse(text = paste0("~0+", arguments[7]))), data=targets))
# Get all combinations of sample groups
group.combinations.mx <- combn(grep(colnames(design), pattern = ind.factors[1], value = T),2)
group.combinations <- foreach(comb = seq(1,ncol(group.combinations.mx)), .combine = c) %do% {paste(group.combinations.mx[,comb], collapse = "-")}
# Rename levels by removing "Sample group and souce prefixes."
colnames(design) <- gsub(colnames(design), pattern = paste0("(", paste(ind.factors, collapse = "|"),")"), replacement = "")
group.combinations <- gsub(group.combinations, pattern = ind.factors[1], replacement = "")
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(contrasts = group.combinations, levels = design)

cat("Creating final models with contrasts...\n")
for(l2ratios in log2ratios.list) {
# fit the linear model 
fit <- lmFit(get(l2ratios), design)
# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
# Look at the numbers of DM CpGs at FDR < 0.05

sink(paste("Analysis summary", suffix, "txt", sep = "."), append = T)
cat("\nNumbers of CpG sites differentially methylated (adjust.method = BH, adj.p.value = 0.05, lfc = 0) between the analyzed groups", if(l2ratios == "log2ratios") {"(original methylation data):"} else {"(binarized methylation data):"},"\n\n")
print(summary(decideTests(fit2))) # adjust.method = "BH", p.value = 0.05, lfc = 0
sink()

# Get a table of differentially methylated probes.
DMPs <- foreach(coef = seq(1, ncol(fit2))) %dopar% {topTable(fit2, num=Inf, coef = coef, genelist = EPICannoSub, adjust.method = "BH", p.value = 0.05, sort.by = "p")}
names(DMPs) <- colnames(fit2)

# Add hg38 coordinates for each CpG site
f.hg38_coords <- function(DMP) {
  gr1 <- makeGRangesFromDataFrame(df = DMP, start.field = "pos", end.field = "pos")
  hg38.coordinates <- as.data.frame(liftOver(gr1, chain = ch19_38))[,c("group_name","start"), drop = F]
  merged.df <- merge(x = DMP, y = hg38.coordinates, by.x = 0, by.y = "group_name", all = T, sort = F)
  colnames(merged.df)[1] <- "Name"
  colnames(merged.df)[3] <- "pos_hg19"
  colnames(merged.df)[54] <- "pos_hg38"
  rownames(merged.df) <- merged.df[["Name"]]
  merged.df[,c(1:3,54,4,6:53)]
}

DMPs <- lapply(DMPs, f.hg38_coords)
assign(paste("DMPs.list", l2ratios, sep = "."), value = DMPs)

# Save the results to a xls workbook.
wb <- paste("wb", l2ratios, sep = ".")
assign(wb, value = createWorkbook())

if(l2ratios == "log2ratios") {datatype <- "original_data"} else {datatype <- "binarized_data"}
for(sheet in names(DMPs)) {
  sheet.name <- substr(paste(sheet, "DMPs", datatype, sep = "."), 1,31)
  tryCatch(expr = {
    addWorksheet(wb = get(wb), sheetName = sheet.name)
  }, error = function(e) {
    sheet.name <<- stri_reverse(substr(stri_reverse(paste(sheet, "DMPs", datatype, sep = ".")),1,31))
    addWorksheet(wb = get(wb), sheetName = sheet.name)
  }
  )
  writeData(x = DMPs[[sheet]], wb = get(wb), sheet = sheet.name, rowNames = F)
}

# Plot the top 4 most significantly differentially methylated CpGs 
if(l2ratios == "log2ratios"){
pdf.name <- paste("4 most significantly differentially methylated CpGs before binarization", suffix, "pdf", sep = ".")} else {
pdf.name <- paste("4 most significantly differentially methylated CpGs after binarization", suffix, "pdf", sep = ".")}

pdf(file = pdf.name, title = pdf.name)
for(group in names(DMPs)) {
  par(mfrow=c(2,2))
  sapply(rownames(DMPs[[group]])[1:4], function(cpg){
    plotCpg(get(l2ratios), cpg=cpg, pheno=targets[[ind.factors[1]]], ylab = "M-values")
  })
  par(mfrow=c(1,1))
  mtext(paste("Group:", sub(group, pattern = "-", replacement = " vs ")), side = 3, line = -1, outer = T)
}
dev.off()
}

# Upset diagrams for DMPs

cat("Drawing upset diagrams for DMPs...\n")
pdf.name <- paste("Upset_plots_for_DMPs", suffix, "pdf", sep = ".")
pdf(title = pdf.name, file = pdf.name, width = 13, height = 7)
for(l2ratios in log2ratios.list) {
  DMPs <- get(paste("DMPs.list", l2ratios, sep = "."))
  if(l2ratios == "log2ratios") {datatype <- "original data set"} else {datatype <- "binarized data set"}
  print(upset(fromList(lapply(DMPs, rownames)), nsets = length(DMPs), nintersects = NA, order.by = c("freq"), decreasing = c(TRUE,FALSE), main.bar.color = cols25()[2], sets.bar.color = cols25()[1]))
  grid.text(paste("Upset plot for DMPs, ", datatype),x = 0.65, y=0.95, gp=gpar(fontsize=12))
}  
dev.off()

# Gene ontology analysis
cat("Performing gene ontology analysis...\n")
for(l2ratios in log2ratios.list) {
  DMPs <- get(paste("DMPs.list", l2ratios, sep = "."))
  if(l2ratios == "log2ratios") {datatype <- "original_data"} else {datatype <- "binarized_data"}
  for(name in names(DMPs)) {
    # Get the significant CpG sites at less than 5% FDR
    sigCpGs <- as.character(DMPs[[name]]$Name[DMPs[[name]]$adj.P.Val<0.05])
    # Get all the CpG sites used in the analysis to form the background
    all <- rownames(mSetSqFlt)
    gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, collection = "GO", plot.bias=F, array.type = "EPIC", sig.genes = F)
    topGO.df <- topGSA(gst, number = 100)
    assign(paste("GO", l2ratios, name, sep = "."), value = topGO.df)
    wb <- paste("wb", l2ratios, sep = ".")
    sheet.name <- substr(paste("GO", l2ratios, name, sep = "."), 1,31)
    tryCatch(expr = {
    addWorksheet(wb = get(wb), sheetName = sheet.name)
    }, error = function(e) {
      sheet.name <<- stri_reverse(substr(stri_reverse(paste("GO", l2ratios, name, sep = ".")),1,31))
      addWorksheet(wb = get(wb), sheetName = sheet.name)
      }
    )
    writeData(x = get(paste("GO", l2ratios, name, sep = ".")), wb = get(wb), sheet = sheet.name)
  }
}

for(l2ratios in log2ratios.list) {
  # Print 100 best CpG plots
  DMPs <- get(paste("DMPs.list", l2ratios, sep = "."))
  if(l2ratios == "log2ratios") {datatype <- "original_data"} else {datatype <- "binarized_data"}
  
  cat("Drawing 100 best CpG plots for ", datatype, "...\n", sep = "")
  foreach(name = names(DMPs)) %dopar% {
    sigCpGs.100 <- DMPs[[name]]$Name[DMPs[[name]]$adj.P.Val<0.05][1:100]
    sigCpGs.100 <- sigCpGs.100[!is.na(sigCpGs.100)]
    DMPs.df <- DMPs[[name]][sigCpGs.100,c("Name", "UCSC_RefGene_Name","adj.P.Val", "logFC")]
    DMPs.df[,2] <- unlist(lapply(str_split(DMPs.df[,2], pattern = ";"), function(x){x[1]}))
    
    if(l2ratios == "log2ratios") {Beta.vals <- Beta.values} else {Beta.vals <- Beta.values.bin}
    pdf.name <- paste("100_most_significant_CpGs", name, datatype, suffix, "pdf", sep = ".")
    if(ncol(Beta.vals) < 63) {width <- 7} else {width <- ceiling(ncol(Beta.vals)/9)}
    pdf(title = pdf.name, file = pdf.name, width = width)
    for(cpgsite in DMPs.df$Name) {
      if(!all(colnames(Beta.vals) == targets$ID)) {stop("The colnames and target IDs do not match.")}
      if(nchar(DMPs.df[DMPs.df$Name == cpgsite,2]) == 0) {gene.name <- "unknown"} else {gene.name <- unlist(strsplit(DMPs.df[DMPs.df$Name == cpgsite,2], split = ";"))[1]}
      padj <- formatC(DMPs.df[DMPs.df$Name == cpgsite,"adj.P.Val"], digits = 3)
      FC <- round(2**(DMPs.df[DMPs.df$Name == cpgsite,"logFC"]),3)
      df1 <- as.data.frame(cbind(Beta.vals[cpgsite,], targets[[ind.factors[1]]]))
      colnames(df1) <- c("Beta-values", ind.factors[1])
      df1[["Beta-values"]] <- as.numeric(df1$`Beta-values`)
      
      plot1 <- ggplot(df1) + 
        geom_point(mapping = aes(x = reorder(targets[[snames]], as.numeric(factor(targets[[ind.factors[1]]]))), y = `Beta-values`, color = targets[[ind.factors[1]]])) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(title = paste("Methylation status of the", cpgsite, "CpG site in gene:", gene.name), x = "Sample names", color = ind.factors[1], subtitle = paste0(name, ": BH-adjusted p-value = ", padj, "; FC = ", FC)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.subtitle = element_text(hjust = 0.5)) +
        scale_color_manual(values=as.vector(cols25()))
      print(plot1)
    }
    dev.off()
  }
  all.sig.cgs <- unique(foreach(name = names(DMPs), .combine = c) %do% {rownames(DMPs[[name]])})
  sink(file = "/dev/null")
  sig.cpgs.mx <- foreach(name = names(DMPs), .combine = cbind) %:%
    foreach(sig.cg = all.sig.cgs, .combine = c) %dopar% {
      if(sig.cg %in% DMPs[[name]]$Name) {
        if(DMPs[[name]][sig.cg,"logFC"] <0) {"Down"} else
        if(DMPs[[name]][sig.cg,"logFC"] >0) {"Up"}} else {
      "Non_sig"}
    }
  sink()
  colnames(sig.cpgs.mx) <- names(DMPs)
  rownames(sig.cpgs.mx) <- all.sig.cgs
  
  sig.cpgs.mx <- cbind(cpgData.annotated[rownames(sig.cpgs.mx),], sig.cpgs.mx)
  sig.cpgs.mx <- sig.cpgs.mx[order(rownames(sig.cpgs.mx)),]

  wb <- paste("wb", l2ratios, sep = ".")
  sheet.name <- substr(paste("Sig_CpGs_summary", datatype, sep = "."), 1,31)
  tryCatch(expr = {
    addWorksheet(wb = get(wb), sheetName = sheet.name)
  }, error = function(e) {
    sheet.name <<- stri_reverse(substr(stri_reverse(paste("Sig_CpGs_summary", datatype, sep = ".")),1,31))
    addWorksheet(wb = get(wb), sheetName = sheet.name)
  }
  )
  writeData(x = sig.cpgs.mx, wb = get(wb), sheet = sheet.name, rowNames = T)
  assign(sheet.name, value = sig.cpgs.mx, envir = .GlobalEnv)

# Order annotation data and Beta values.
if(l2ratios == "log2ratios") {Beta.vals <- Beta.values} else {Beta.vals <- Beta.values.bin}
if(!all(rownames(Beta.vals) %in% rownames(cpgData.annotated))) {stop("Not all CpG sites have been found in Illumina annotation data.")}
cpgData.annotated <- cpgData.annotated[match(rownames(Beta.vals), rownames(cpgData.annotated)),]
Beta.vals <- Beta.vals[, match(targets$ID, colnames(Beta.vals))]

# Performing the CpGs methylation analysis in the context of genes, regions, and strands.
cat("Performing the CpGs methylation analysis in the context of genes, regions, and strands (", datatype, ")...\n", sep = "")

level.names <- levels(anno[[ind.factors[1]]])
comb.list <- strsplit(group.combinations, split = "-")
names(comb.list) <- group.combinations

if(length(grep(unique.genes, pattern = "\\(b\\)")) >0) {stop("Some genes are assigned to two different strands simultaneously.")}

all.res <- foreach(gene = unique.genes) %dopar% {
  gene.regex <- sub(gene, pattern = "(^.*)\\(([pm])\\)$", replacement = "\\1\\\\(\\2\\\\)")
  strand.type <- sub(gene, pattern = "(^.*)\\(([pm])\\)$", replacement = "\\2")
  if(strand.type == "p") {strand.fin <- "+"} else 
    if(strand.type == "m") {strand.fin <- "-"} else
      {stop("Strand type identification has failed.")}

  cpgData.annotated.filtered <- cpgData.annotated[grepl(cpgData.annotated$Genes, pattern = paste0("(^|;)", gene.regex)) & cpgData.annotated$strand %in% strand.fin, , drop = F]
  if(nrow(cpgData.annotated.filtered) >0) {
    unique.annots <- unique(gsub(subsetByOverlaps(x = subset(annotations, subset = symbol == sub(gene, pattern = "\\([pm]\\)", replacement = "")), ranges = makeGRangesFromDataFrame(cpgData.annotated.filtered))$type, pattern = "(hg19_|genes_|_gencode|_fantom)", replacement = ""))
  all.annots <- foreach(annot = paste0(unique.annots, "(", strand.type, ")")) %do% {
    Beta.vals.subset <- subset(Beta.vals, subset = rownames(Beta.vals) %in% rownames(cpgData.annotated.filtered) & grepl(cpgData.annotated$Annotations, pattern = sub(annot, pattern = "\\(([pm])\\)", replacement = "\\\\(\\1\\\\)")))
    cpgs.no <- nrow(Beta.vals.subset)
    cpgsNames <- paste(rownames(Beta.vals.subset), collapse = "&")
    if(cpgs.no >0) {
    merged.data <- merge(x = as.data.frame(colMeans(Beta.vals.subset)), y = anno, by.x = 0, by.y = 0, sort = F)
    colnames(merged.data)[2] <- "colMeans"
    means.SDs <- unlist(with(merged.data, by(merged.data, INDICES = get(ind.factors[1]), FUN = function(x) {return(c(mean(x[["colMeans"]]), sd(x[["colMeans"]])))})))
    
    kwallis.res <- with(merged.data, kruskal.test(colMeans ~ get(ind.factors[1])))
      wilcox.res <- foreach(comb = names(comb.list)) %do% {
        merged.data.comb <- subset(merged.data, subset = get(ind.factors[1]) %in% comb.list[[comb]])
        with(merged.data.comb, wilcox.test(colMeans ~ get(ind.factors[1])))
      }
      final.res <- c(cpgs.no, cpgsNames, means.SDs, kwallis.res$p.value, sapply(wilcox.res, FUN = function(x) {x$p.value}))
    } else {
      final.res <- c(cpgs.no, rep(NA, length(level.names)*2 + length(comb.list)+2))
    }
      names(final.res) <- c("CpGs.no", "CpGs.names",
                            as.character(t(outer(sub(level.names, pattern = "[12]$", replacement = ""), c("average.betas.mean", "average.betas.sd"), FUN = paste, sep = "."))),
                            "all.groups.Kruskal-Wallis.test.pval",
                            paste(names(comb.list), "Wilcoxon.rank.sum.test.pval", sep = "."))
      final.res
  }
  names(all.annots) <- paste0(unique.annots, "(", strand.type, ")")
  all.annots
  }
}
names(all.res) <- unique.genes
all.res <- all.res[!sapply(all.res, is.null)]
all.res.names <- foreach(index = seq(all.res), .combine = c) %do% {paste(names(all.res)[index], names(all.res[[index]]), sep = ".")}

all.res.df <- as.data.frame(t(as.data.frame(all.res, make.names = NA, check.names = F)))
rownames(all.res.df) <- all.res.names
all.res.df <- all.res.df[all.res.df$CpGs.no >0, , drop = F]
genes.regions <- data.frame(gene.name = sub(rownames(all.res.df), pattern = "(^.*)\\.([^\\.]*$)", replacement = "\\1"), gene.region = sub(rownames(all.res.df), pattern = "(^.*)\\.([^\\.]*$)", replacement = "\\2"))
all.res.df <- all.res.df %>% tibble::add_column(genes.regions, .before = 1)
all.res.df <- all.res.df[order(as.numeric(all.res.df[["all.groups.Kruskal-Wallis.test.pval"]])),]
kwallis.mins <- cbind(with(all.res.df, by(all.res.df, INDICES = gene.name, function(x) {min(as.numeric(x[["all.groups.Kruskal-Wallis.test.pval"]]))})))
kwallis.mins <- kwallis.mins[match(unique(all.res.df$gene.name), rownames(kwallis.mins)),]
kwallis.mins <- kwallis.mins[kwallis.mins <0.05]
all.res.df <- all.res.df[all.res.df$gene.name %in% names(kwallis.mins), , drop = F]
final.order <- foreach(name = names(kwallis.mins), .combine = c) %dopar% {which(all.res.df$gene.name %in% name)}
all.res.df <- all.res.df[final.order,]
all.res.df <- all.res.df %>% tibble::add_column(index = seq(1, nrow(all.res.df)), .before = 1)
addWorksheet(wb = get(wb), sheetName = "CpGs.GeneRegions.analysis")
writeData(x = all.res.df, wb = get(wb), sheet = "CpGs.GeneRegions.analysis", rowNames = F)

genes.list <- unique(all.res.df$gene.name)
invisible(foreach(gene.no = seq(length(genes.list))) %dopar% {
  gene <- genes.list[gene.no]
  strand.type <- sub(gene, pattern = "(^.*)\\(([pm])\\)$", replacement = "\\2")
  if(strand.type == "p") {strand.fin <- "+"} else 
    if(strand.type == "m") {strand.fin <- "-"} else
    {stop("Strand type identification has failed.")}
  
  pdf.name <- paste("GeneRegions.analysis", sprintf(fmt = "%07d", gene.no), datatype, "pdf.tmp", sep = ".")
  pdf(title = pdf.name, file = pdf.name, width = 10)
  all.res.df.sub <- all.res.df[all.res.df$gene.name == gene, , drop = F]
  
  invisible(foreach(region = all.res.df.sub$gene.region) %do% {
    all.res.df.sub.sub <- all.res.df.sub[all.res.df.sub$gene.region == region, , drop = F]
    Beta.vals.subs <- subset(Beta.vals, subset = grepl(cpgData.annotated$Genes, pattern = paste0("(^|;)", sub(gene, pattern = "\\(([pm])\\)", replacement = "\\\\(\\1\\\\)"), "(;|$)")) & 
                               grepl(cpgData.annotated$Annotations, pattern = paste0("(^|;)", sub(region, pattern = "\\(([pm])\\)", replacement = "\\\\(\\1\\\\)"), "(;|$)")) &
                               cpgData.annotated$strand == strand.fin)
    if(nrow(all.res.df.sub.sub) != 1 | all.res.df.sub.sub$CpGs.no != nrow(Beta.vals.subs)) {stop("Beta values subsetting by the gene, region and strand failed.")}
    Aver.Betas.df <- data.frame(Average.Betas = colMeans(Beta.vals.subs))
    Aver.Betas.df.merged <- merge(x = Aver.Betas.df, y = anno, by.x = 0, by.y = 0)
    all.res.df.sub.sub <- setNames(as.data.frame(t(all.res.df.sub.sub[grep(colnames(all.res.df.sub.sub), pattern = "\\.pval$")])), nm = "p-values")
    all.res.df.sub.sub$`p-values` <- formatC(as.numeric(all.res.df.sub.sub$`p-values`, digits = 3))
    all.res.df.sub.sub <- tibble::add_column(.data = all.res.df.sub.sub, tests = rownames(all.res.df.sub.sub), .before = 1)
    all.res.df.sub.sub$tests <- sub(gsub(all.res.df.sub.sub$tests, pattern = "\\.", replacement = " "), pattern = " pval$", replacement = "")
    
    box.plot <- ggplot(Aver.Betas.df.merged) +
      aes(x = get(ind.factors[1]), y = Average.Betas, fill = get(ind.factors[1])) +
      geom_violin() +
      geom_boxplot(width = 0.1) +
      geom_jitter(color = "azure4", size = 1) +
      stat_summary(geom = "point", fun = mean, fill = "yellow", shape = 24) +
      annotate(geom = "table", x = (length(level.names)+1)/2, y = max(Aver.Betas.df$Average.Betas) * 1.3, label = list(all.res.df.sub.sub), table.theme = ttheme_gtminimal, vjust = 1, hjust = 0.5) +
      scale_fill_manual(values = as.character(cols25())) +
      labs(title = paste("Comparison of beta values distribution, gene:", gene, ", region:", region),
           x = ind.factors[1], y = paste("Average beta values")) +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
    print(box.plot)
  })
  dev.off()
})

pdf.files <- list.files(pattern = paste0("^GeneRegions\\.analysis\\.[0-9]{7}\\.", datatype, "\\.pdf\\.tmp$"))
pdf_combine(pdf.files, output = paste("GeneRegions.analysis", datatype, suffix, "pdf", sep = "."))
unlink(pdf.files)

# Merge betas with annotations.
cpgData.annotated.meth <- GRanges(cpgData.annotated,
                                  Annotations = cpgData.annotated$Annotations,
                                  Genes = cpgData.annotated$Genes,
                                  betas = Beta.vals)

# Differential methylation analysis of regions
DMRs.Fisher.list <- list()
ranges.mean.beta.vals.list <- list()

for(name in names(DMPs)) {
cat("Performing differential methylation analysis of regions for ", datatype, ", comparison: ", name, "...\n", sep = "")  

myAnnotation <- cpg.annotate(object = get(l2ratios), datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = name, arraytype = "EPIC")

if(!all(names(myAnnotation@ranges) == rownames(cpgData.annotated))) {stop("Rownames between CpG annotations and CpG data do not match.")}
myAnnotation@ranges@strand <- Rle(cpgData.annotated$strand)

DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
DMRs.name <- paste("DMRs", l2ratios, name, sep = ".")
assign(DMRs.name, value = DMRs)

results.ranges <- extractRanges(dmrcoutput = DMRs, genome = genome.ver)
rr.name <- paste("results.ranges", l2ratios, name, sep = ".")
assign(rr.name, value = results.ranges)

results.ranges.sig <- results.ranges[results.ranges$Fisher<0.05]
rr.name.sig <- paste("results.ranges.sig", l2ratios, name, sep = ".")
assign(rr.name.sig, value = results.ranges.sig)

# Calculate Fisher statistics for stranded DMR data.
cat("Performing the DMR analysis for individual strands (", datatype, "), comparison: ", name, "...\n", sep = "")

dmrcate.Fisher.stranded <- function (object, lambda = 1000, C = NULL, pcutoff = "fdr", consec = FALSE, 
                              conseclambda = 10, betacutoff = NULL, min.cpgs = 2) 
{
  stopifnot(is(object, "CpGannotated"))
  stopifnot(lambda >= 1)
  stopifnot(pcutoff == "fdr" | (0 <= pcutoff & pcutoff <= 1))
  stopifnot(C >= 0.2)
  if (consec & is.null(conseclambda)) {
    stop("Consecutive CpG bandwidth must be specified")
  }
  object <- data.frame(ID = names(object@ranges), strand = strand(object@ranges), weights = abs(object@ranges$stat), 
                       CHR = seqnames(object@ranges), pos = start(object@ranges), 
                       diff = object@ranges$diff, indfdr = object@ranges$ind.fdr, 
                       is.sig = object@ranges$is.sig)
  object <- object[order(object$CHR, object$pos), ]
  if (is.null(C) & !consec) {
    C = 2
  }
  if (consec) {
    lambda = conseclambda
    message(paste("Consecutive mode specified, lambda is now set at", 
                  conseclambda, "consecutive CpGs."))
    if (is.null(C)) {
      stop("Error: argument C must be specified (in CpG sites) for consecutive mode.")
    }
    object$realcoordforconsec <- object$pos
    object$pos <- unlist(sapply(as.numeric(table(object$CHR)), 
                                function(x) 1:x))
  }
  lag = lambda
  chr.unique <- unique(c(as.character(object$CHR)))
  fitted <- lapply(chr.unique, fitParallel, object = object, 
                   consec = consec, conseclambda = conseclambda, lambda = lambda, 
                   C = C)
  object <- rbind.fill(fitted)
  object$fdr <- p.adjust(object$raw, method = "BH")
  if (pcutoff == "fdr") {
    nsig <- sum(object$is.sig)
    if (nsig == 0) {
      txt <- "The FDR you specified in cpg.annotate() returned no significant CpGs, hence there are no DMRs.\n    Try specifying a value of 'pcutoff' in dmrcate() and/or increasing 'fdr' in cpg.annotate()."
      stop(paste(strwrap(txt, exdent = 2), collapse = "\n"))
    }
    pcutoff <- sort(object$fdr)[nsig]
  }
  object$sig <- (object$fdr <= pcutoff)
  if (nrow(object) == 0) {
    txt <- "No signficant regions found. Try increasing the value of\n    'pcutoff' in dmrcate() and/or 'fdr' in cpg.annotate()."
    stop(paste(strwrap(txt, exdent = 2), collapse = "\n"))
  }
  message("Demarcating regions for individual strands...")
  chr.N <- as.character(object$CHR)
  pos.N <- object$pos
    sig.N <- object$sig
    N <- length(sig.N)
    n.K <- which(sig.N)
    K <- length(n.K)
    stopifnot(K >= 2)
    pos.K <- pos.N[n.K]
    chr.K <- chr.N[n.K]
    jump_chr.k <- (chr.K[-1] != chr.K[-K])
    jump_pos.k <- (diff(pos.K) > lag)
    jump.k <- (jump_chr.k | jump_pos.k)
    ksegments.A2 <- Segment(jump.k)
    A <- nrow(ksegments.A2)
    kstart.A <- ksegments.A2[, "start"]
    kend.A <- ksegments.A2[, "end"]
    realpos.K <- pos.K
    if (consec) {
      realpos.N <- object$realcoordforconsec
      realpos.K <- realpos.N[n.K]
    }
    start.A <- realpos.K[kstart.A]
    end.A <- realpos.K[kend.A]
    chr.A <- chr.K[kstart.A]
    stopifnot(all(chr.K[kend.A] == chr.A))
    fmt <- "%s:%1d-%1d"
    coord.A <- sprintf(fmt, chr.A, start.A, end.A)
    nstart.A <- n.K[kstart.A]
    nend.A <- n.K[kend.A]
    width.A <- nend.A + 1 - nstart.A
    a.Z <- rep(seq(A), width.A)
    fn <- function(a) seq(from = nstart.A[a], to = nend.A[a])
    l.listA <- lapply(seq(A), fn)
    n.Z <- unlist(l.listA)
    region.N <- rep(NA_integer_, N)
    region.N[n.Z] <- a.Z    
    levels <- seq(A)
    region.N <- factor(region.N, levels = levels)
  
  antypes <- c("both", "plus", "minus")
  fisher.res <- foreach(an.type = antypes, .combine = cbind) %dopar% {
    
    REGIONSTAT <- function(field, fn) {
      x.N <- object[[field]]
      x.R <- tapply(x.N, region.N, fn)
      c(x.R)
    }
    fn_Fisher <- function(x) pchisq((sum(log(x)) * -2), df = length(x) * 
                                      2, lower.tail = FALSE)
    
    if (an.type == "both") {
      strand.set = c("+", "-")
    } else
      if (an.type == "plus") {
        strand.set = c("+")
      } else
        if(an.type == "minus") {
          strand.set = c("-")
        }
    region.N[!object$strand %in% strand.set] <- NA_integer_
    dmr.chr <-c(tapply(X = chr.N, INDEX = region.N, FUN = unique))
    dmr.ranges <-c(tapply(X = pos.N, INDEX = region.N, FUN = function(x) {paste(range(x), collapse = "-")}))
    dmr.coords <- paste(dmr.chr, dmr.ranges, sep = ":")
#    dmr.strands <-c(tapply(X = object$strand, INDEX = region.N, FUN = function(x) {paste(unique(x), collapse = "&")}))
    dmr_cpgs.no <- c(table(region.N)) 
    dmr.stranded.res <- data.frame(coord = dmr.coords, no.CpGs = dmr_cpgs.no, Fisher = REGIONSTAT("indfdr", fn_Fisher))
    setNames(dmr.stranded.res, nm = paste(colnames(dmr.stranded.res), an.type, sep = "."))
    }
  fisher.res <- fisher.res[fisher.res$no.CpGs.both >= min.cpgs,]
  fisher.res <- fisher.res[order(fisher.res$Fisher.both, -fisher.res$no.CpGs.both),]
  rownames(fisher.res) <- NULL
  
  message("Done!")
  return(new("data.frame", fisher.res))
}

DMRs.Fisher <- dmrcate.Fisher.stranded(myAnnotation, lambda=1000, C=2)

if(nrow(DMRs.Fisher) != length(results.ranges)) {stop("The numbers of DMRs in the stranded and original DMR analysis (comparison: ", name, ") do not match.")}

DMRs.Fisher <- subset(DMRs.Fisher, subset = Fisher.both <0.05 | Fisher.plus <0.05 | Fisher.minus <0.05)
DMRs.Fisher.list <- append(DMRs.Fisher.list, list(DMRs.Fisher))

# DNAseI hypersensitive sites
dnase <- read.table(file = dhs.bed, sep = "\t", header = F, stringsAsFactors = F)
dnaseData <- GRanges(seqnames=dnase[,1],
                     ranges=IRanges(start=dnase[,2], end=dnase[,3]),
                     strand=Rle(rep("*",nrow(dnase))),
                     data=dnase[,5])

# TFBS
tfbs <- read.table(file = tfbs.bed, sep = "\t", header = F, stringsAsFactors = F)
tfbsData <- GRanges(seqnames=tfbs[,1],
                     ranges=IRanges(start=tfbs[,2], end=tfbs[,3]),
                     strand=Rle(rep("*",nrow(tfbs))),
                     data=tfbs[,5])

# Draw all significant DMRs in the genomic context 
if(length(get(rr.name.sig)) > 0) {
cat("Drawing all significant DMRs in the genomic context for ", datatype, ", comparison: ", name, "...\n", sep = "")
  if(threads >25){threads <- 25
  registerDoMC(threads)}
invisible(foreach(range = seq(1, length(get(rr.name.sig)))) %dopar% {
chrom <- as.character(seqnames(get(rr.name.sig)[range]))
start <- as.numeric(start(get(rr.name.sig)[range]))
end <- as.numeric(end(get(rr.name.sig)[range]))
# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))

iTrack <- IdeogramTrack(genome = genome.ver, chromosome = chrom)
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
# rTrack <- UcscTrack(genome=genome.ver, chromosome=chrom, track="NCBI RefSeq", 
#                     from=minbase, to=maxbase, trackType="GeneRegionTrack", 
#                     rstarts="exonStarts", rends="exonEnds", gene="name", 
#                     symbol="name2", transcript="name", strand="strand", 
#                     fill="lightgreen",stacking="squish", name="RefSeq", 
#                     showId=TRUE, geneSymbol=TRUE)

# Extract data on CpGs in DMR
GencodeV10.dmr <- subsetByOverlaps(GencodeV10, get(rr.name.sig)[range])
cpgData.dmr <- subsetByOverlaps(cpgData.annotated.meth, get(rr.name.sig)[range])
GencodeV10.dmr$gene_id <- GencodeV10.dmr$gene_name
txdb.dmr <- makeTxDbFromGRanges(GencodeV10.dmr)

# Remove excessive data from DMR
strand(cpgData.dmr) <- "*"
cpgData.dmr$Annotations <- NULL
cpgData.dmr$Genes <- NULL

rTrack <- GeneRegionTrack(range = txdb.dmr, 
                           start = minbase, 
                           end = maxbase, 
                           name = "RefSeq\nhg19",
                           fill = "lightgreen",
                           stacking = "squish",
                          col.title="black", col.axis="black",
                          fontcolor.item = "black",
                          cex.item = 0.6,
                          rotation.item = 0)

# methylation data track
methTrack <- DataTrack(range=cpgData.dmr, groups=targets[[ind.factors[1]]], genome = genome.ver,
                       chromosome=chrom, ylim=c(-0.05,1.05), col=cols25(),
                       type=c("a","p"), name="DNA Meth.\n(beta value)",
                       background.panel="white", legend=TRUE,
                       col.title="black", col.axis="black")
# DNaseI hypersensitive site data track
dnaseTrack <- DataTrack(range=dnaseData, genome=genome.ver, name="DHSS", 
                        type="gradient", chromosome=chrom, ylim = range(dnase[,5]),
                        col.title="black", col.axis="black")

# TFBS
tfbsTrack <- DataTrack(range=tfbsData, genome=genome.ver, name="TFBS", 
                        type="gradient", chromosome=chrom, ylim = range(tfbs[,5]),
                       gradient = c("white","red"),
                       col.title="black", col.axis="black")

# DMR position data track
dmrTrack <- AnnotationTrack(start=start, end=end, genome=genome.ver, name="DMR", 
                            chromosome=chrom,fill="darkred",
                            col.title="black")
tracks <- list(iTrack, gTrack, dmrTrack, methTrack, dnaseTrack, tfbsTrack,
               rTrack)
sizes <- c(1,2,1,5,2,2,4) # set up the relative sizes of the tracks

pdf.name <- paste("DMR_context_plot", range, name, datatype, "pdf.tmp", sep = ".")
tryCatch(expr = {
pdf(title = pdf.name, file = pdf.name)
plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks), 
           main = paste0("DMR ", range, "\nGenes: ", paste0(unique(GencodeV10.dmr$gene_id), collapse = ", ")),
           cex.main = 1, shape = "arrow",
           exonAnnotation = "gene")},
error = function(e) {
pdf(title = pdf.name, file = pdf.name, height = 14)
  plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
             add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks), 
             main = paste0("DMR ", range, "\nGenes: ", paste0(unique(GencodeV10.dmr$gene_id), collapse = ", ")),
             cex.main = 1, shape = "arrow",
             exonAnnotation = "gene")})
dev.off()
})
threads <- as.numeric(arguments[5])
registerDoMC(threads)
}
if(length(get(rr.name.sig))>0) {
pdffiles <- foreach(range = seq(1, length(get(rr.name.sig))), .combine = c) %do% {
  paste("DMR_context_plot", range, name, datatype, "pdf.tmp", sep = ".")}
pdf_combine(input = pdffiles, output = paste("DMR_context_plots", name, datatype, suffix, "pdf", sep = "."))
unlink(pdffiles)}

cat("Calculating mean beta values within every DMR for each sample (", datatype, "), comparison: ", name, "...\n", sep = "")
ranges.mean.beta.vals <- foreach(range = seq(1, length(get(rr.name))), .combine = rbind) %dopar% {
    rr.subset <- as.data.frame(subsetByOverlaps(cpgData.annotated.meth, get(rr.name)[range]))
    rr.subset.plus <- rr.subset[rr.subset$strand == "+", , drop = F]
    rr.subset.minus <- rr.subset[rr.subset$strand == "-", , drop = F]
    mx.tmp <- as.matrix(rbind(
    cbind(rr.subset %>% dplyr::summarize(across(.cols = c(1:7), .fns = function(x) {paste(unique(unlist(strsplit(as.character(x), split = ";"))), collapse = "&")})),
          rr.subset %>% dplyr::summarize(across(.cols = !c(1:7), .fns = function(x) {mean((x))}))),
    if(nrow(rr.subset.plus) >0 & nrow(rr.subset.minus) >0) {
    rbind(
    cbind(rr.subset.plus %>% dplyr::summarize(across(.cols = c(1:7), .fns = function(x) {paste(unique(unlist(strsplit(as.character(x), split = ";"))), collapse = "&")})),
          rr.subset.plus %>% dplyr::summarize(across(.cols = !c(1:7), .fns = function(x) {mean((x))}))),
    cbind(rr.subset.minus %>% dplyr::summarize(across(.cols = c(1:7), .fns = function(x) {paste(unique(unlist(strsplit(as.character(x), split = ";"))), collapse = "&")})),
          rr.subset.minus %>% dplyr::summarize(across(.cols = !c(1:7), .fns = function(x) {mean((x))})))
    )
      }
    ))
    rownames(mx.tmp) <- rep(paste0("DMR_", range), nrow(mx.tmp))
    mx.tmp
    }
colnames(ranges.mean.beta.vals) <- sub(colnames(ranges.mean.beta.vals), 
                                           pattern = "betas\\.", 
                                           replacement = "")
rownames(ranges.mean.beta.vals) <- foreach(row = seq(1, nrow(ranges.mean.beta.vals)), .combine = c)  %do% {paste(rownames(ranges.mean.beta.vals)[row], ranges.mean.beta.vals[row, 5], sep = "_")}
ranges.mean.beta.vals[,"start"] <- sapply(strsplit(ranges.mean.beta.vals[,"start"], split = "&"), FUN = function(y) {paste(range(y), collapse = "-")})
colnames(ranges.mean.beta.vals)[2:3] <- c("range", "CpGs.positions")

sheet.name <- substr(paste("DMR_mean.betas", name, datatype, sep = "."), 1,31)
tryCatch(expr = {
  addWorksheet(wb = get(wb), sheetName = sheet.name)
}, error = function(e) {
  sheet.name <<- stri_reverse(substr(stri_reverse(paste("DMR_mean.betas", name, datatype, sep = ".")),1,31))
  addWorksheet(wb = get(wb), sheetName = sheet.name)
}
)
writeData(x = ranges.mean.beta.vals, wb = get(wb), sheet = sheet.name, rowNames = T)
ranges.mean.beta.vals.list <- append(ranges.mean.beta.vals.list, list(t(ranges.mean.beta.vals)))
}
names(DMRs.Fisher.list) <- names(DMPs)
names(ranges.mean.beta.vals.list) <- names(DMPs)

assign(x = paste("DMRs.Fisher.list", l2ratios, sep = "."), value = DMRs.Fisher.list)
assign(x = paste("ranges.mean.beta.vals.list", l2ratios, sep = "."), value = ranges.mean.beta.vals.list)

cat("Performing common DMRs identification and analysis for ", datatype, "...\n", sep = "")

vec.union <- sort(unique(unlist(sapply(DMRs.Fisher.list, function(x) {return(c(paste(x[["coord.both"]],"both", sep = "#"), paste(x[["coord.plus"]], "plus", sep = "#"), paste(x[["coord.minus"]], "minus", sep = "#")))}))))
vec.union <- sub(vec.union, pattern = ":", replacement = "#")

ranges.mean.beta.vals.list <- sapply(ranges.mean.beta.vals.list, FUN = function(x) {x["strand",] <- sub(sub(sub(x["strand",], pattern = "^\\+$", replacement = "plus"), 
                                                                              pattern = "^\\-$", replacement = "minus"), 
                                                                          pattern = "^(\\-\\&\\+|\\+\\&\\-)$", replacement = "both"); return(x)})

ranges.mean.beta.vals.list.vec <- foreach(df.no = seq(1, length(ranges.mean.beta.vals.list))) %:%
  foreach(i = seq(1, ncol(ranges.mean.beta.vals.list[[df.no]])), .combine = c) %do% {
    paste(as.character(ranges.mean.beta.vals.list[[df.no]][c(1,2,5,6,7), i]), collapse = "#")}
names(ranges.mean.beta.vals.list.vec) <- names(DMPs)

ranges.mean.beta.vals.list.vec.unique <- unique(unlist(ranges.mean.beta.vals.list.vec))

vec.union.full <- mclapply(vec.union, FUN = function(x) {res <- grep(ranges.mean.beta.vals.list.vec.unique, 
                                                              pattern = x, value = T)
                                                              if(length(res) > 1) {
                                                                stop("The unequivocal identification of full names for some unique DMRs has failed.")
                                                              } else {return(res)}}, mc.cores = threads)
if(any(sapply(vec.union.full, class) == "try-error")) {stop("The unequivocal identification of full names for some unique DMRs has failed.")}
vec.union.full <- unlist(vec.union.full)

common.DMRs <- foreach(test.val = seq(1, length(vec.union.full)), .combine = rbind) %dopar% {
  bool.vec <- foreach(test.vec = names(DMPs), .combine = c) %do% {
    if(!all(rownames(ranges.mean.beta.vals.list[[test.vec]])[-c(1:7)] %in% rownames(anno))) {stop("The sample names do not match annotations.")}
    vec.union.full[test.val] %in% ranges.mean.beta.vals.list.vec[[test.vec]]
  }
  if(sum(bool.vec) == 0) {stop("The identification of the correct DMR group with statistically significant methylation changes has failed.")}
  
  min.vec <- min(which(bool.vec == 1))
  
  if(ncol(ranges.mean.beta.vals.list[[min.vec]][, as.data.frame(t(ranges.mean.beta.vals.list[[min.vec]][c(1,2,5,6,7),])) %>% 
                                                tidyr::unite("united", 1:5, sep = "#") == vec.union.full[test.val], drop = F]) != 1) 
  {stop(paste("Some ranges in the", names(ranges.mean.beta.vals.list)[min.vec], "comparison are missing or are duplicated."))}
  
  merged.df <- merge(x = ranges.mean.beta.vals.list[[min.vec]][, as.data.frame(t(ranges.mean.beta.vals.list[[min.vec]][c(1,2,5,6,7),])) %>% 
                                                                 tidyr::unite("united", 1:5, sep = "#") == vec.union.full[test.val], drop = F][-c(1:7), , drop = F], y = anno, by.x = 0, by.y = 0)
  merged.df[,2] <- as.numeric(merged.df[,2])
  merged.df <- merged.df %>% dplyr::select(2, ind.factors[1])
  mean.vec <- unlist(as.list(with(merged.df, by(merged.df[,1], INDICES = get(ind.factors[1]), FUN = function(x) {c(mean(x),sd(x))}))))
  final.df <- rbind(c(bool.vec, mean.vec))
  rownames(final.df) <- vec.union.full[test.val]
  final.df
}
if(!all(unique(sub(colnames(common.DMRs)[(length(DMPs)+1):(ncol(common.DMRs))], pattern = "[12]$", replacement = "")) == levels(anno[[ind.factors[1]]]))) 
  {stop("Column names in the common.DMRs table and independent factor levels in the anno table do not match.")}
colnames(common.DMRs) <- c(paste("DMRcate", names(DMPs), sep = "."), as.character(t(outer(levels(anno[[ind.factors[1]]]), c(".Mean_meth_freq.", ".SD_meth_freq."), FUN = paste0))))
common.DMRs <- as.data.frame(common.DMRs)
common.DMRs <- common.DMRs[order(rownames(common.DMRs)),]

unique.coords.stranded <- lapply(strsplit(rownames(common.DMRs), split = "#"), FUN = function(x) {return(c(sprintf("%s:%s", x[1], x[2]), x[3]))})

Fisher.res.calc <- function(bool) {bool <- as.logical(bool)
    foreach(unique.coord.no = seq(unique.coords.stranded), .combine = rbind) %dopar% {
  unique.coords <- unique.coords.stranded[[unique.coord.no]]
  names(unique.coords) <- paste(c("coord", "Fisher"), unique.coords[2], sep = ".")
  foreach(DMRs.Fisher.name = DMRs.Fisher.list, .combine = c) %do% {
    Fisher.test.res <- subset(DMRs.Fisher.name, 
           subset = get(names(unique.coords)[1]) == unique.coords[1])[,names(unique.coords)[2]]
    if(bool) {
    if(length(Fisher.test.res) == 1) {if(Fisher.test.res < 0.05) {"TRUE"} else {"FALSE"}} else
              if(length(Fisher.test.res) == 0) {"FALSE"} else 
              {stop("Identification of the corresponding Fisher result has failed.")}} else {
                if(length(Fisher.test.res) == 1) {Fisher.test.res} else
                  if(length(Fisher.test.res) == 0) {NA} else 
                  {stop("Identification of the corresponding Fisher result has failed.")}
                }
  }
}
}

common.DMRs.Fisher.full <- Fisher.res.calc(bool = F)
colnames(common.DMRs.Fisher.full) <- paste("Fisher_padj", names(DMPs), sep = ".")
common.DMRs.Fisher.bool <- as.character(!(is.na(common.DMRs.Fisher.full) | common.DMRs.Fisher.full >= 0.05))
common.DMRs[,1:length(DMPs)] <- common.DMRs.Fisher.bool
common.DMRs <- common.DMRs %>% add_column(as.data.frame(common.DMRs.Fisher.full), .after = length(DMPs))
common.DMRs <- common.DMRs[order(rowSums(sapply(common.DMRs[,1:length(DMPs)], as.logical)), decreasing = T),]
common.DMRs <- common.DMRs %>% tibble::add_column(Fisher_padj.min.val = dplyr::select(., paste("Fisher_padj", names(DMPs), sep = ".")) %>% as.matrix() %>% rowMins(na.rm = T), .after = length(DMPs)*2)

sheet.name <- substr(paste("Common.DMRs", datatype, sep = "."), 1,31)
tryCatch(expr = {
  addWorksheet(wb = get(wb), sheetName = sheet.name)
}, error = function(e) {
  sheet.name <<- stri_reverse(substr(stri_reverse(paste("Common.DMRs", datatype, sep = ".")),1,31))
  addWorksheet(wb = get(wb), sheetName = sheet.name)
}
)
writeData(x = common.DMRs, wb = get(wb), sheet = sheet.name, rowNames = T, keepNA = T, na.string = "NA")
assign(paste("common.DMRs", l2ratios, sep = "."), value = common.DMRs)

common.DMRs <- common.DMRs[common.DMRs[,1:length(DMPs), drop = F] %>% sapply(., FUN = as.logical) %>% rowSums(.) >0, , drop = F]
if(max(common.DMRs$Fisher_padj.min.val) >= 0.05) {stop("Adjusted p-values' filtering has failed.")}
common.DMRs <- common.DMRs[order(common.DMRs$Fisher_padj.min.val),]

common.DMRs.maxes <- common.DMRs %>% dplyr::transmute(across(.cols = 1:length(DMPs), .fns = as.logical)) %>% t() %>% as.data.frame() %>% sapply(., which.max)

common.DMRs.mean.betas <- foreach(DMR = rownames(common.DMRs), .combine = cbind) %dopar% {
  coords <- unlist(strsplit(DMR, split = "#"))
  comp <- common.DMRs.maxes[names(common.DMRs.maxes) == DMR]
  matching.col <- ranges.mean.beta.vals.list[[comp]][,(ranges.mean.beta.vals.list[[comp]]["seqnames",] == coords[1]) & (ranges.mean.beta.vals.list[[comp]]["range",] == coords[2]) & (ranges.mean.beta.vals.list[[comp]]["strand",] == coords[3]), drop = F]
  if(ncol(matching.col) != 1) {stop("The identification of ", DMR, " has failed.")} else
  {
  matching.col <- t(as.data.frame(t(matching.col)) %>% dplyr::transmute(across(.cols = -c(1:7), .fns = as.numeric)))
  if(nrow(matching.col) != nrow(anno)) {stop("The number of rows in the ", DMR, " column does not match the number of annotated samples.")}
  colnames(matching.col) <- DMR
  matching.col }
}

if(!all(rownames(common.DMRs.mean.betas) == rownames(anno))) {stop("Sample names do not match annotations for some tested DMRs.")}
if(!all(colnames(common.DMRs.mean.betas) == rownames(common.DMRs))) {stop("The order/number of unique DMR names does not match.")}

rownames(common.DMRs.mean.betas) <- sub(rownames(common.DMRs.mean.betas), pattern = paste0("^(", paste(levels(anno[[ind.factors[1]]]), collapse = "|"), ")\\."), replacement = "")
common.DMRs.mean.betas <- common.DMRs.mean.betas[order(rownames(common.DMRs.mean.betas)),]
common.DMRs.mean.betas <- common.DMRs.mean.betas %>% as.data.frame() %>% add_column(Sample.names = rownames(.), .before = 1)

write.table(x = common.DMRs.mean.betas, row.names = F, sep = ";", file = paste("Unique_DMRs_mean.beta.vals", datatype, suffix, "csv", sep = "."))

sink(paste("Analysis summary", suffix, "txt", sep = "."), append = T)
cat("\nNumbers of differentially methylated DMRs (both strands, adjust.method = BH, adj.p.value = 0.05, lfc = 0) between the analyzed groups", if(l2ratios == "log2ratios") {"(original methylation data):"} else {"(binarized methylation data):"},"\n\n")

rr.list <- ls(pattern = paste0("^results\\.ranges\\.sig\\.", l2ratios))
f.sig.dmrs <- function(x) {v.tmp <- as.vector(table(get(x)$meandiff >0, useNA = "no"))
v.tmp[3] <- sum(v.tmp)
names(v.tmp) <- c("Down", "Up", "Total")
v.tmp}

sig.dmrs.df <- sapply(rr.list, f.sig.dmrs)
colnames(sig.dmrs.df) <- sub(colnames(sig.dmrs.df), pattern = paste0("^results\\.ranges\\.sig\\.", l2ratios, "\\."), replacement = "")
print(sig.dmrs.df)

sink()

saveWorkbook(wb = get(wb), overwrite = T, file = paste("Methylation analysis results", datatype, suffix, "xlsx", sep = "."))
}

# save(list = c("mSetRaw", "mSetSq", "mSetSqFlt", "rgSet", "detP", "detP.tibble"), file = paste("Methylation analysis preliminary data", suffix, "RData", sep = "."))
rm(mSetRaw, mSetSq, mSetSqFlt, rgSet, detP, detP.tibble)

save.image <- function(file){save(list=grep(ls(all.names = TRUE, envir = .GlobalEnv), pattern = "^arguments$", value = T, invert = T), file = file)}

save.image(paste("Methylation analysis results", suffix, "RData", sep = "."))
} else
{
  cat("Loading the existing RData object...\n")
  load(paste("Methylation analysis results", suffix, "RData", sep = "."))
  read.args()
}
  if(any(sel.CpGs != "NA")) {
    cat("Drawing plots visualizing frequency of DNA methylation for selected CpG sites...\n")
    
    for(l2ratios in log2ratios.list) {
    if(l2ratios == "log2ratios") {
      DMPs <- DMPs.list.log2ratios
      datatype <- "original_data" } else
    {
      DMPs <- DMPs.list.log2ratios.bin
      datatype <- "binarized_data"
    }
    
    if(l2ratios == "log2ratios") {Beta.vals <- Beta.values} else {
      Beta.vals <- Beta.values.bin}
      
    pdf.name <- paste(basename(CpGs.signature.file), "listed_CpG_sites", datatype, suffix, "pdf", sep = ".")
    if(ncol(Beta.vals) < 63) {width <- 7} else {width <- ceiling(ncol(Beta.vals)/9)}
    pdf(title = pdf.name, file = pdf.name, width = width)
    
    for(CpG in sel.CpGs) {
      labels <- mapply(function(DMPs, name, CpG) {df.tmp <- DMPs[DMPs$Name == CpG,c("logFC", "adj.P.Val")]
      if(nrow(df.tmp) == 1) {
        paste0(name, ": BH−adjusted p−value = ", formatC(df.tmp[,"adj.P.Val"], digits = 3), "; FC = ", round(2**(df.tmp[,"logFC"]),3))} else
      if(nrow(df.tmp) == 0){
        paste0(name, ": insignificant result")} else {
      stop(paste("The analyzed CpG site", CpG, "occured more than once in the result table."))}}, DMPs = DMPs, name = names(DMPs), CpG = CpG)
      labels <- gsub(labels, pattern = "−", replacement = "-")
      
    sel.CpG <- cpgData.annotated[rownames(cpgData.annotated) %in% CpG,]
      
    if(nrow(sel.CpG) > 1) {stop("CpGs subsetting has failed.")}
    if(nrow(sel.CpG) == 1) {
    
      if(!all(colnames(Beta.vals) == targets$ID)) {stop("The colnames and target IDs do not match.")}
      df1 <- as.data.frame(cbind(Beta.vals[CpG,], targets[[ind.factors[1]]]))
      colnames(df1) <- c("Beta-values", ind.factors[1])
      df1[["Beta-values"]] <- as.numeric(df1[["Beta-values"]])
      
      plot1 <- ggplot(df1) + 
        geom_point(mapping = aes(x = reorder(targets[[snames]], as.numeric(factor(targets[[ind.factors[1]]]))), y = `Beta-values`, color = targets[[ind.factors[1]]])) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(title = paste("Methylation status of the", rownames(sel.CpG), "CpG site, on the", paste0("'", sel.CpG$strand, "'"), "strand,\ngene(s):", sel.CpG$Genes), x = "Sample names", color = ind.factors[1]) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=as.vector(cols25()))
      plot1 <- plot1 + ylim(0,1+(length(labels)*0.05)) + 
        annotate("text", x = (nrow(plot1$data)+1)/2, y = 1+(length(labels)*0.03), label = paste(labels, collapse = "\n"))
      print(plot1)
    } else {
      plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("The", CpG, "site was not found in the 'cpgData.annotated' R object used in the present analysis.\nThis can happen if the CpG name is incorrect or the corresponding CpG site did not pass the preliminary filtering step."))
    }
  } 
    dev.off()
  }
  }

if(!file.exists(paste("Methylation_data", suffix, "soft", sep = "."))) {
  library(wateRmelon)
  library(methylumi)
  all.idats <- list.files(unique(file.path(dataDirectory, targets$Slide)), pattern = ".*\\.idat$", full.names = T)
  idats.keep <- foreach(idat = all.idats, .combine = c) %do% {sjmisc::str_contains(idat, pattern = targets$Basename, logic = "OR")}
  sel.idats <- all.idats[idats.keep]
  barcodes.list <- unique(basename(sub(sel.idats, pattern = "_(Grn|Red)\\.idat$", replacement = "")))
  dataTable <- dataTable[dataTable %>% dplyr::select(all_of(c("Sentrix_ID", "Sentrix_Position"))) %>% unite(col = "United", sep = "_") %>% .[,"United"] %in% barcodes.list,]
  dataTable <- dataTable %>% unite(c("Sentrix_ID", "Sentrix_Position"), col = "barcodes", sep = "_")
  dataTable <- dataTable[order(dataTable$barcodes),]
  dir.create("sel.idats.dir")
  if(!all(file.copy(from = sel.idats, to = "sel.idats.dir", overwrite = T))) {stop("An error occurred while copying idat files to the sel.idats.dir")}
  if(!all(dataTable$barcodes == unique(sub(list.files("sel.idats.dir/"), pattern = "_(Grn|Red)\\.idat$", replacement = "")))) {stop("Barcodes of idat files do not match those privided in the attached data frame.")}

  epic.data <- readEPIC(idatPath = "sel.idats.dir", pdat = dataTable, force = T)
  epic.data.norm <- normalizeMethyLumiSet(epic.data)
  epic.data.MethyLumiM <- as(epic.data, "MethyLumiM")
  epic.data.norm.MethyLumiM <- as(epic.data.norm, "MethyLumiM")
  
  produceGEOSampleInfoTemplate(lumiNormalized = epic.data.norm.MethyLumiM)
  GEOsampleInfo <- fread("GEOsampleInfo.txt")
  unlink("GEOsampleInfo.txt")
  stopifnot(dataTable[["barcodes"]] == GEOsampleInfo[["sampleID"]])
  GEOsampleInfo$Sample_title <- dataTable$Sample_Name
  GEOsampleInfo$Sample_source_name_ch1 <- dataTable$Sample_Source
  GEOsampleInfo$Sample_organism_ch1 <- "Homo sapiens"
  GEOsampleInfo$Sample_characteristics_ch1 <- paste("Group", dataTable$Sample_Group, sep = ": ")
  GEOsampleInfo$Sample_label_ch1 <- "Cy3, Cy5"
  GEOsampleInfo$Sample_description <- dataTable$Sample_Description
  GEOsampleInfo$Sample_platform_id <- "GPL21145"
  GEOsampleInfo$Sample_supplementary_file <- NULL
  sel.idats.names <- list.files("sel.idats.dir/")
  sel.idats.names <- sel.idats.names[foreach(ID = GEOsampleInfo$sampleID, .combine = c) %do% {which(grepl(sel.idats.names, pattern = ID))}]
  supplementary.data <- cbind(sel.idats.names[seq(1,length(sel.idats.names), 2)], sel.idats.names[seq(2,length(sel.idats.names), 2)])
  colnames(supplementary.data) <- c("Sample_supplementary_file", "Sample_supplementary_file")
  GEOsampleInfo <- cbind(GEOsampleInfo, supplementary.data)
  
  
  unlink("sel.idats.dir", recursive = T)
  
}






sessionInfo()
proc.time()
date()
cat("All done.\n")
