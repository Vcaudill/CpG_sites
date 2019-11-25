source("R_scripts/DataFormatting/WTAA_consensus.R")
source("R_scripts/DataFormatting/MUTAA.R")
source("R_scripts/DataFormatting/Drastic_AA_Change.R")
source("R_scripts/DataFormatting/CPG_Function.R")
source("R_scripts/DataFormatting/SynNonSyn.R")

#Cycle through the data folder
Virus_info<- list.files("data/", pattern=".fasta")
for(i in 1:length(Virus_info)){
  print(i)
  source("R_scripts/DataFormatting/Freq.R")
  # must place your file as a txt takes a few minutes 
  viruplace = paste0("data/",Virus_info[i])
  DF<-Freq(viruplace)
  name <- as.character(Virus_info[i])
  truename<-unlist(strsplit(as.character(Virus_info[i]),".fasta"))
 
  DF$wtnt_consensus<-as.character(DF$wtnt_consensus)
  DF$Virus<-(truename)

  DF<-getWTAA(DF)
  DF<-getMUTAA(DF)
  DF<-big_aa_change(DF)
  DF<-CPG_site(DF)
  DF<-synFunction(DF)

  write.csv(DF, file = paste0("output/Csv/",truename,".csv"))
}
