
my.list<- list.files("output/Csv", pattern=".csv")
data_points = data.frame("order"= 1:length(my.list), "Virus"= 0,'total'=0,"AsynNC_C" =0,"AnonsynNC_C"=0,'TsynNC_C'=0, "TnonsynNC_C"= 0,'AsynC_LCLS' = 0, 'AnonsynC_LCLS'=0,'AsynNC_LCLS'= 0,'AnonsynNC_LCLS'=0,'TsynC_LCLS' = 0, 'TnonsynC_LCLS' = 0, 'TsynNC_LCLS' = 0, 'TnonsynNC_LCLS' = 0, 'AsynC_UCLS'= 0, 'AnonsynC_UCLS'= 0, 'AsynNC_UCLS'= 0, 'AnonsynNC_UCLS'= 0, 'TsynC_UCLS'= 0, 'TnonsynC_UCLS'=0, 'TsynNC_UCLS' = 0,'TnonsynNC_UCLS' =0   )
Virus_info<-read.csv("output/Extra_CpG_info.csv")
final_t<-data.frame("File"=character(),"t_test_a"=integer(),"t_test_t"=integer(),"ancestor_t_test_a"=integer(),"ancestor_t_test_t"=integer())
Virus_info<-Virus_info[!(names(Virus_info) %in% c("X" ))]
for (i in 1:length(my.list)){  

data<-read.csv(paste0("output/Csv/",my.list[i]))
tt<-Virus_info[grep(strsplit(my.list[i],".csv"), Virus_info$File), ]

tt$t_test_a<-t.test(data$Freq[data$wtnt_consensus =="a" & data$makesCpG == 1],mu=0, alternative = "greater")$p.value
tt$t_test_t<-t.test(data$Freq[data$wtnt_consensus =="t" & data$makesCpG == 1],mu=0, alternative = "greater")$p.value


tt$ancestor_t_test_a<-t.test(data$aFreq[data$ancestor =="a" & data$ancestor_makesCpG == 1],mu=0, alternative = "greater")$p.value
tt$ancestor_t_test_t<-t.test(data$aFreq[data$ancestor =="t" & data$ancestor_makesCpG == 1],mu=0, alternative = "greater")$p.value

final_t <- rbind(final_t, tt[c("File","t_test_a","t_test_t","ancestor_t_test_a","ancestor_t_test_t")])


if(my.list[length(my.list)]==my.list[i]){
   Virus_info<- merge(Virus_info, final_t, by= "File")
   write.csv(Virus_info, file = "output/Extra_CpG_info.csv")
}

}


boxplot(data$aFreq[data$ancestor =="a" & data$ancestor_makesCpG == 1])
boxplot(data$aFreq[data$ancestor =="t" & data$ancestor_makesCpG == 1])
