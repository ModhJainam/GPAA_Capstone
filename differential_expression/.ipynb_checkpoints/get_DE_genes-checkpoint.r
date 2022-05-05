library('DESeq2')
library(devtools)
library(ashr)
#no adjustments
gcountData_AA <- as.matrix(read.csv("gene_count_matrix_AA.csv", row.names='gene_id', check.names=FALSE))
colData_AA <- read.csv('GPAA_AA_b13_metadata_fixed.csv', row.names = 1)

reorder_idx <- match(colnames(gcountData_AA), rownames(colData_AA))
colData_AA <- colData_AA[reorder_idx, ]

colData_AA$binary_pathology <- factor(colData_AA$binary_pathology)

noconfound_gdds <-DESeqDataSetFromMatrix(countData = gcountData_AA, colData = colData_AA, design = ~binary_pathology)
noconfound_gdds <- estimateSizeFactors(noconfound_gdds)
noconfound_gdds <- DESeq(noconfound_gdds)

noconfound_gres2<-results(noconfound_gdds, contrast=c("binary_pathology","1","0"))
noconfound_gres2_ashr <-lfcShrink(dds=noconfound_gdds, contrast=c('binary_pathology', "1", "0"), type="ashr")
CAD_noconfound_selected_AA <- rownames(subset(noconfound_gres2_ashr,padj<=0.05 & abs(log2FoldChange)>=1))

########################################################################################################

gcountData_LAD <- as.matrix(read.csv("gene_count_matrix_LAD.csv", row.names='gene_id', check.names=FALSE))
colData_LAD <- read.csv('GPAA_LAD_b13_metadata_fixed.csv', row.names = 1)

reorder_idx <- match(colnames(gcountData_LAD), rownames(colData_LAD))
colData_LAD <- colData_LAD[reorder_idx, ]

colData_LAD$binary_pathology <- factor(colData_LAD$binary_pathology)

noconfound_gdds <-DESeqDataSetFromMatrix(countData = gcountData_LAD, colData = colData_LAD, design = ~binary_pathology)
noconfound_gdds <- estimateSizeFactors(noconfound_gdds)
noconfound_gdds <- DESeq(noconfound_gdds)

noconfound_gres2<-results(noconfound_gdds, contrast=c("binary_pathology","1","0"))
noconfound_gres2_ashr <-lfcShrink(dds=noconfound_gdds, contrast=c('binary_pathology', "1", "0"), type="ashr")
CAD_noconfound_selected_LAD <- rownames(subset(noconfound_gres2_ashr,padj<=0.05 & abs(log2FoldChange)>=1))

########################################################################################################
#adjust for confounding
gcountData_AA <- as.matrix(read.csv("gene_count_matrix_AA.csv", row.names='gene_id', check.names=FALSE))
colData_AA <- read.csv('GPAA_AA_b13_metadata_fixed.csv', row.names = 1)

reorder_idx <- match(colnames(gcountData_AA), rownames(colData_AA))
colData_AA <- colData_AA[reorder_idx, ]

colData_AA$binary_pathology <- factor(colData_AA$binary_pathology)

asr_gdds <-DESeqDataSetFromMatrix(countData = gcountData_AA, colData = colData_AA, design = ~age+white+male+binary_pathology)
asr_gdds <- estimateSizeFactors(asr_gdds)
asr_gdds <- DESeq(asr_gdds)

asr_gres2<-results(asr_gdds, contrast=c("binary_pathology","1","0"))
asr_gres2_ashr <-lfcShrink(dds=asr_gdds, contrast=c('binary_pathology', "1", "0"), type="ashr")
CAD_asr_selected_AA <- rownames(subset(asr_gres2_ashr,padj<=0.05 & abs(log2FoldChange)>=1))

########################################################################################################

gcountData_LAD <- as.matrix(read.csv("gene_count_matrix_LAD.csv", row.names='gene_id', check.names=FALSE))
colData_LAD <- read.csv('GPAA_LAD_b13_metadata_fixed.csv', row.names = 1)

reorder_idx <- match(colnames(gcountData_LAD), rownames(colData_LAD))
colData_LAD <- colData_LAD[reorder_idx, ]

colData_LAD$binary_pathology <- factor(colData_LAD$binary_pathology)

asr_gdds <-DESeqDataSetFromMatrix(countData = gcountData_LAD, colData = colData_LAD, design = ~age+white+male+binary_pathology)
asr_gdds <- estimateSizeFactors(asr_gdds)
asr_gdds <- DESeq(asr_gdds)

asr_gres2<-results(noconfound_gdds, contrast=c("binary_pathology","1","0"))
asr_gres2_ashr <-lfcShrink(dds=asr_gdds, contrast=c('binary_pathology', "1", "0"), type="ashr")
CAD_asr_selected_LAD <- rownames(subset(asr_gres2_ashr,padj<=0.05 & abs(log2FoldChange)>=1))

########################################################################################################

writeLines(CAD_noconfound_selected_AA, con = 'noconfound_AA_DE_genes.txt', sep="\n")
writeLines(CAD_noconfound_selected_LAD, con='noconfound_LAD_DE_genes.txt', sep="\n")

writeLines(CAD_asr_selected_AA, con = 'asr_AA_DE_genes.txt', sep="\n")
writeLines(CAD_asr_selected_LAD, con='asr_LAD_DE_genes.txt', sep="\n")