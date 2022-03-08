library(dplyr)
library(caret)
library(dprep)
library(readr)
library(stringr)

AA_genomic <- read_csv("/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/preprocessing/GPAA_AA_genomic_b13.csv")
LAD_genomic <- read_csv("/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/preprocessing/GPAA_LAD_genomic_b13.csv")

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

AA_genomic.lessnzv.no_sampleid.lesscor.norm <- znorm(AA_genomic.lessnzv.no_sampleid.lesscor)
LAD_genomic.lessnzv.no_sampleid.lesscor.norm <- znorm(LAD_genomic.lessnzv.no_sampleid.lesscor)

rownames(AA_genomic.lessnzv.no_sampleid.lesscor.norm) = AA_genomic$sample_id
rownames(LAD_genomic.lessnzv.no_sampleid.lesscor.norm) = LAD_genomic$sample_id

write.csv(AA_genomic.lessnzv.no_sampleid.lesscor.norm, "GPAA_AA_genomic_b13_preprocessed.csv", row.names=TRUE)
write.csv(LAD_genomic.lessnzv.no_sampleid.lesscor.norm, "GPAA_LAD_genomic_b13_preprocessed.csv", row.names=TRUE)

##############################################################################################

AA_genomic.lessnzv.no_sampleid.norm <- znorm(AA_genomic.lessnzv.no_sampleid)
LAD_genomic.lessnzv.no_sampleid.norm <- znorm(LAD_genomic.lessnzv.no_sampleid)

rownames(AA_genomic.lessnzv.no_sampleid.norm) = AA_genomic$sample_id
rownames(LAD_genomic.lessnzv.no_sampleid.norm) = LAD_genomic$sample_id

asr_AA_DE_genes = readLines('/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/differential_expression/asr_AA_DE_genes.txt')
asr_LAD_DE_genes = readLines('/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/differential_expression/asr_LAD_DE_genes.txt')
noconfound_AA_DE_genes = readLines('/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/differential_expression/noconfound_AA_DE_genes.txt')
noconfound_LAD_DE_genes = readLines('/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/differential_expression/noconfound_LAD_DE_genes.txt')

asr_AA_DE_genes<-str_extract(string=asr_AA_DE_genes, pattern="[^|]+")
asr_LAD_DE_genes<-str_extract(string=asr_LAD_DE_genes, pattern="[^|]+")
noconfound_AA_DE_genes<-str_extract(string=noconfound_AA_DE_genes, pattern="[^|]+")
noconfound_LAD_DE_genes<-str_extract(string=noconfound_LAD_DE_genes, pattern="[^|]+")

asr_AA_towrite = AA_genomic.lessnzv.no_sampleid.norm[names(AA_genomic.lessnzv.no_sampleid.norm)[names(AA_genomic.lessnzv.no_sampleid.norm) %in% asr_AA_DE_genes]]
asr_LAD_towrite = LAD_genomic.lessnzv.no_sampleid.norm[names(LAD_genomic.lessnzv.no_sampleid.norm)[names(LAD_genomic.lessnzv.no_sampleid.norm) %in% asr_LAD_DE_genes]]
noconfound_AA_towrite = AA_genomic.lessnzv.no_sampleid.norm[names(AA_genomic.lessnzv.no_sampleid.norm)[names(AA_genomic.lessnzv.no_sampleid.norm) %in% noconfound_AA_DE_genes]]
noconfound_LAD_towrite = LAD_genomic.lessnzv.no_sampleid.norm[names(LAD_genomic.lessnzv.no_sampleid.norm)[names(LAD_genomic.lessnzv.no_sampleid.norm) %in% noconfound_LAD_DE_genes]]

write.csv(asr_AA_towrite, "GPAA_AA_genomic_b13_preprocessed_asr_DE_genes.csv", row.names=TRUE)
write.csv(asr_LAD_towrite, "GPAA_LAD_genomic_b13_preprocessed_asr_DE_genes.csv", row.names=TRUE)
write.csv(noconfound_AA_towrite, "GPAA_AA_genomic_b13_preprocessed_noconfound_DE_genes.csv", row.names=TRUE)
write.csv(noconfound_LAD_towrite, "GPAA_LAD_genomic_b13_preprocessed_noconfound_DE_genes.csv", row.names=TRUE)