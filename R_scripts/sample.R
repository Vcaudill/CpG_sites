library(ape)
library(seqinr) 
#install.packages("phylotools")
library(phylotools) 
library(dplyr)

viruplace<-"data/EnterovirusA_VP2_EVA71.fasta"
name<-"EnterovirusA_VP2_EVA71"
DataSet <- read.fasta(viruplace)


viruplace<-"data/HBV_A_Polymerase.fasta"
name<-"HBV_A_Polymerase"
DataSet <- read.fasta(viruplace)

viruplace<- "data/HBV_A_S.fasta"
name<- "HBV_A_S"
DataSet <- read.fasta(viruplace)

viruplace<-"data/HBV_B_Polymerase.fasta"
name<-"HBV_B_Polymerase.fasta"
DataSet <- read.fasta(viruplace)

viruplace<-"data/HBV_B_PreC-Core.fasta"
name<-"HBV_B_PreC-Core.fasta"
DataSet <- read.fasta(viruplace)

viruplace<-"data/HBV_C_polymerase.fasta"
name<-"HBV_C_polymerase.fasta"
DataSet <- read.fasta(viruplace)

viruplace<-"data/HBV_C_PreC-Core.fasta"
name<-"HBV_C_PreC-Core.fasta"
DataSet <- read.fasta(viruplace)

viruplace<-"data/HBV_C_S.fasta"
name<-"HBV_C_S.fasta"
DataSet <- read.fasta(viruplace)

viruplace<-"data/Parainfulenza3_CDS.fasta"
name<-"Parainfulenza3_CDS"
DataSet <- read.fasta(viruplace)


sample_size= 200 #sequences
DataSet <- read.fasta(viruplace)
# DF<-as.data.frame(DataSet)
nrow(DataSet)
if (nrow(DataSet) > sample_size) {
    DS <- dplyr::sample_n(DataSet, size = sample_size, replace=F)
    size = sample_size
}else{
    DS <- dplyr::sample_n(DataSet, size = as.numeric(nrow(DataSet)))
}

sampleviruplace = paste("/Users/victoria/Desktop/",name,"_sample.fasta",sep="")
dat2fasta(DS, outfile = sampleviruplace)


