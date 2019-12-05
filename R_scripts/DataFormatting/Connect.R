
source("R_scripts/DataFormatting/Freq.R")
source("R_scripts/DataFormatting/WTAA_consensus.R")
source("R_scripts/DataFormatting/MUTAA.R")
source("R_scripts/DataFormatting/Drastic_AA_Change.R")
source("R_scripts/DataFormatting/CPG_Function.R")
source("R_scripts/DataFormatting/SynNonSyn.R")

ancestor<-read.csv("output/Final_CpG_List.csv")
ancestor2<-data.frame("File"=character(),"ca"=integer(),"tg"=integer(),"ancestor_ca"=integer(),"ancestor_tg"=integer())
#Cycle through the data folder
Virus_info<- list.files("data/", pattern=".fasta")

for(i in 1:length(Virus_info)){
  print(Virus_info[i])
  ancestor_line<-ancestor[grep(Virus_info[i], ancestor$File), ]
  viruplace = paste0("data/",Virus_info[i])
  DF<-Freq(viruplace, ancestor_line$ancestor_seq)
  name <- as.character(Virus_info[i])
  truename<-unlist(strsplit(as.character(Virus_info[i]),".fasta"))
 
  DF$wtnt_consensus<-as.character(DF$wtnt_consensus)
  DF$ancestor<-as.character(DF$ancestor)
  DF$Virus<-(truename)

  DF<-getWTAA(DF)
  DF<-getMUTAA(DF)
  DF<-big_aa_change(DF)
  DF <-CPG_site(DF)
  DF<-synFunction(DF)
  
  ancestor_line$ca<- DF$ca[1]
  ancestor_line$tg<- DF$tg[1]
  ancestor_line$ancestor_ca<- DF$ancestor_ca[1]
  ancestor_line$ancestor_tg<- DF$ancestor_tg[1]
  ancestor2 <- rbind(ancestor2, ancestor_line[c("File","ca","tg","ancestor_ca","ancestor_tg")])
  
  DF<-DF[!(names(DF) %in% c("ca","tg","ancestor_ca", "ancestor_tg" ))]
  write.csv(DF, file = paste0("output/Csv/",truename,".csv"))
  if(Virus_info[length(Virus_info)]==Virus_info[i]){
      ancestor<- merge(ancestor, ancestor2, by= "File")
      write.csv(ancestor, file = "output/Extra_CpG_info.csv")
  }
 
}
