install.packages(c("dplyr", "caret", "readr", "stringr"))
library(dplyr)
library(caret)
#library(dprep)
library(readr)
library(stringr)

AA_genomic <- read_csv("/project/gpaa/machine_learning/jainam_capstone/GPAA_samples_AA_batches0-15_FPKM.csv")
LAD_genomic <- read_csv("/project/gpaa/machine_learning/jainam_capstone/GPAA_samples_LAD_batches0-15_FPKM.csv")

nzv_AA <- nearZeroVar(AA_genomic)
nzv_LAD <- nearZeroVar(LAD_genomic)
AA_genomic.lessnzv <- AA_genomic[,-nzv_AA]
LAD_genomic.lessnzv <- LAD_genomic[,-nzv_LAD]

AA_genomic.lessnzv.no_sampleid <- subset(AA_genomic.lessnzv, select = -c(sample_id))
LAD_genomic.lessnzv.no_sampleid <- subset(LAD_genomic.lessnzv, select = -c(sample_id))

correlations_AA <- cor(AA_genomic.lessnzv.no_sampleid)
correlations_LAD <- cor(LAD_genomic.lessnzv.no_sampleid)
highcorrelations_AA <- findCorrelation(correlations_AA, cutoff = .80)
highcorrelations_LAD <- findCorrelation(correlations_LAD, cutoff = .80)
AA_genomic.lessnzv.no_sampleid.lesscor <- AA_genomic.lessnzv.no_sampleid[,-highcorrelations_AA]
LAD_genomic.lessnzv.no_sampleid.lesscor <- LAD_genomic.lessnzv.no_sampleid[,-highcorrelations_LAD]

#AA_genomic.lessnzv.no_sampleid.lesscor.norm <- znorm(AA_genomic.lessnzv.no_sampleid.lesscor)
#LAD_genomic.lessnzv.no_sampleid.lesscor.norm <- znorm(LAD_genomic.lessnzv.no_sampleid.lesscor)

AA_genomic.lessnzv.no_sampleid.lesscor.norm <- scale(AA_genomic.lessnzv.no_sampleid.lesscor)
LAD_genomic.lessnzv.no_sampleid.lesscor.norm <- scale(LAD_genomic.lessnzv.no_sampleid.lesscor)

rownames(AA_genomic.lessnzv.no_sampleid.lesscor.norm) = AA_genomic$sample_id
rownames(LAD_genomic.lessnzv.no_sampleid.lesscor.norm) = LAD_genomic$sample_id

write.csv(AA_genomic.lessnzv.no_sampleid.lesscor.norm, "GPAA_AA_samples_preprocessed.csv", row.names=TRUE)
write.csv(LAD_genomic.lessnzv.no_sampleid.lesscor.norm, "GPAA_LAD_samples_preprocessed.csv", row.names=TRUE)

##############################################################################################
