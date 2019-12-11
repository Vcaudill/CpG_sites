
#Cycle through the csv folder
Virus_info<- list.files("output/Csv", pattern=".csv")
Virus_info_names<-read.csv("output/Final_CpG_List.csv")

source("R_scripts/Tables/MakeWilcoxTables_ancestor.R")
source("R_scripts/graphs/M_frequency_ancestor_graph.R")

for(i in 1:length(Virus_info)){
 
  Name<-Virus_info_names[grep(strsplit(Virus_info[i],".csv"), Virus_info_names$File), ]
  truename<-unlist(strsplit(as.character(Virus_info[i]),".csv"))
  table_output<- "output/supplenentary_based_on_ancestor/WilcoxTables/"
  graph_output<-"output/supplenentary_based_on_ancestor/M_frequency_graphs/"
  #table
  DF=Tables(truename, paste0('output/Csv/',Virus_info[i]))
  Pvalues=Wilcox_test(DF, truename)#get error x must be numeric
  makeTable(Pvalues, truename,  Name$name, table_output)
  #M freq graph
  comparing_CpG_Syn_Nonsyn_new(truename, Name$name, paste0('output/Csv/',Virus_info[i]), graph_output)
  #dev.off()
  
}
