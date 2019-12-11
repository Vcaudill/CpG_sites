# Calculate dinucleotide frequencies
library(msa)
library(ape)
library(car)
library(seqinr)
library(purrr)
library(tidyverse)

#dir.create("output/Dinuc/")
Fastafiles<-list.files("data/",pattern=".fasta$")

for (i in 1:length(Fastafiles)){
    fas1<-read.alignment(paste0("data/",Fastafiles[i]),format = "fasta")
    
    freq<-list()
    Rho<-list()
    Zscores<-list()
    for (j in 1: length(fas1[[3]])){
        seq<-fas1[[3]][[j]]
        seq<-substring(seq, 1:nchar(seq),1:nchar(seq))
        
        #record dinucleotide freq
        freq[[j]]<-as.data.frame(seqinr::count(seq,2))
        names(freq)[j]<-fas1[[2]][j]
        Rho[[j]]<-as.data.frame(rho(seq,2))
        names(Rho)[j]<-fas1[[2]][j]
        Zscores[[j]]<-as.data.frame(zscore(seq, modele="base"))
        names(Zscores)[j]<-fas1[[2]][j]
    }
    
    
    dint<-data.frame(diNT=freq[[1]][1])
    Rho.values<-data.frame(diNT=freq[[1]][1])
    Z.scores<-data.frame(diNT=freq[[1]][1])
    for (j in 1:length(freq)) {
        dt1<-freq[[j]]
        dt2<-Rho[[j]]
        dt3<-Zscores[[j]]
        dint<-cbind(dint,dt1[,2])
        Rho.values<-cbind(Rho.values, dt2[,2])
        Z.scores<-cbind(Z.scores,dt3[,2])
    }
    
    
    
    tb1<-data.frame(diNT=dint[,1])
    tb1$Freq<-rowMeans(dint[,2:ncol(dint)])
    tb1$Rho<-rowMeans(Rho.values[,2:ncol(Rho.values)], na.rm=1)
    tb1$Z<-rowMeans(Z.scores[,2:ncol(Z.scores)], na.rm=1)
    
    write.csv(tb1,paste0("output/consensus/Dinuc/Dinuc.freq.",substr(Fastafiles[i], 1,nchar(Fastafiles[i])-6),".csv"))
    
}




