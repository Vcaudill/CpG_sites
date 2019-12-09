library(seqinr) 
ancestor<-read.csv("output/Final_CpG_List.csv")
ancestor2<-data.frame("File"=character(),"ca"=integer(),"tg"=integer(),"ancestor_ca"=integer(),"ancestor_tg"=integer())
#Cycle through the data folder
Virus_info<- list.files("output/Csv/", pattern=".csv")
Virus_fasta<- list.files("data/", pattern=".fasta")
for(i in 1:length(Virus_info)){
    print(Virus_info[i])
    print(Virus_fasta[i])
    virus_basic <- read.fasta(paste0("data/",Virus_fasta[i]))
    number_of_seqs <- length(virus_basic)
    virus_csv<-read.csv(paste0("output/Csv/",Virus_info[i]))
    num_of_nuc<-nrow(virus_csv)
    ancestor$Number_of_Sequences[i]<-number_of_seqs
    ancestor$Number_of_Nucleotides[i]<-num_of_nuc
}

write.csv(ancestor, file = "output/Final_CpG_List.csv")
