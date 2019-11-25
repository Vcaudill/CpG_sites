# your floder where all pages are saved
#setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")


#loop using a csv to find the file name, strating and stoping points, and readign frame

Virus_info<- read.csv("data/list/Final_CpG_list.csv")
# DataSet <-read.fasta("DengueVirus1.fasta_pruned.mu.trim05.txt")
for(i in 1:nrow(Virus_info)){
  print(i)
  source("R_scripts/Basic_Bioinfo_Rscripts/Freq.R")
  # must place your file as a txt takes a few minutes 
  #setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/virus")
  viruplace = paste(Virus_info$Fasta_File_Path[i],"/fasta/", Virus_info$File[i], sep="")
  DF<-Freq(viruplace)
  name = as.character(Virus_info$File[i])
  splitname<-unlist(strsplit(as.character(Virus_info$File[i]),".fasta"))
  truename<-splitname[1]
  DF$wtnt_consensus<-as.character(DF$wtnt_consensus)
  DF$Virus<-(truename)
  #setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")
  
  #  
  # make sure that M is the first amino acid check all virus
  #no longer needed
  #start = Virus_info$start[i]
  #end = Virus_info$stop[i]
  #DF<-DF[c(start:(nrow(DF)-end)),]
  #DF$num<-(1:nrow(DF))
  
  
  source("R_scripts/Basic_Bioinfo_Rscripts/WTAA_consensus.R")
  DF<-getWTAA(DF)
  
  source("R_scripts/Basic_Bioinfo_Rscripts/MUTAA.R")
  DF<-getMUTAA(DF)
  
  source("R_scripts/Basic_Bioinfo_Rscripts/Drastic_AA_Change.R")
  DF<-big_aa_change(DF)
  
  source("R_scripts/Basic_Bioinfo_Rscripts/CPG_Function.R")
  DF<-CPG_site(DF)
  
  source("R_scripts/Basic_Bioinfo_Rscripts/SynNonSyn.R")
  DF<-synFunction(DF)
  
  #How to save data
  #truenameRda= paste('data/Rda/', truename, ".Rda", sep="")
  #save("virus" = DF,file= truenameRda)
  truenameCSV= paste(Virus_info$Fasta_File_Path[i],"/Csv/", truename,".csv", sep="")
  write.csv(DF, file = truenameCSV)
}

 source("R_scripts/Basic_Bioinfo_Rscripts/SynNonSyn.R")
DF<-synFunction(DF)

source("R_scripts/Basic_Bioinfo_Rscripts/CPG_Function.R")
DF<-CPG_site(DF)
#How to load data
virusname = 'EnterovirusB_VP2.fasta_pruned.mu.trim05.txt'
splitname<-unlist(strsplit(virusname,".fasta"))
truename<-splitname[1]

#truenameRda= paste('data/Rda/', truename, ".Rda", sep="")
truenameCSV= paste('data/data_2019/data_used/Csv/', truename, ".csv", sep="")
#write.csv(DF, file = "MyData.csv")

#setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-/Rda_Files")
# load(truenameRda)
###############################################
# graphs/ tables 
#setwd("~/Desktop/Git/CpG/Something_Cool-CpG_Sites-")
# source("R_scripts/BioInfo/Synonymous and nonsynonymous graph.R")
# SynNonsynAT(DF)
# SynNonsynCG(DF)
# 
# 
# source("R_scripts/BioInfo/MakesCPG Graph.R")
# CPGNoCPGAT(DF)
# 
# source("R_scripts/BioInfo/CPG_Syn_Nonsyn_graph.R")
# comparing_CpG_Syn_Nonsyn (DF)

# Wilcox test
#source("stats.R")
#Wilcox_test(DF)
# that tables
#source("R_scripts/BioInfo/RyanWilcox.R")
#Wilcox_test(DF, truename)

#source("R_scripts/graphs/redoplot.R")
#comparing_CpG_Syn_Nonsyn_new(DF, truename)

#########
########


