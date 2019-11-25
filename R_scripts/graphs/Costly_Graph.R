
library(scales)
library(plotrix)
library(sfsmisc)
#add viruses with names and csv roots
#update my.list and name.list by alpha order and DNA or RNA
#abline (ln230)
#fill virus info into excel sheet
#_______
BKpolyomavirus_VP1<-read.csv("data/data_2019/data_used/Csv/BKpolyomavirus_VP1.CSV")

DengueVirus1 <- read.csv("data/data_2019/data_used/Csv/DengueVirus1.CSV")
DengueVirus2<- read.csv("data/data_2019/data_used/Csv/DengueVirus2.CSV")
DengueVirus3<-read.csv("data/data_2019/data_used/Csv/DengueVirus3.CSV")
DengueVirus4<-read.csv("data/data_2019/data_used/Csv/DengueVirus4.CSV")
EnterovirusA_VP1<-read.csv("data/data_2019/data_used/Csv/EnterovirusA_VP1_EVA71.CSV")
EnterovirusA_VP2<-read.csv("data/data_2019/data_used/Csv/EnterovirusA_VP2_EVA71.CSV")
EnterovirusB_VP1<-read.csv("data/data_2019/data_used/Csv/EnterovirusB_VP1_Echovirus30.CSV")
EnterovirusB_VP2<-read.csv("data/data_2019/data_used/Csv/EnterovirusB_VP2.CSV")
EnterovirusC_VP1<-read.csv("data/data_2019/data_used/Csv/EnterovirusC_VP1_Polio2.CSV")
EnterovirusC_VP2<-read.csv("data/data_2019/data_used/Csv/Enterovirus_C_VP2_Polio2.CSV")
EnterovirusD_VP1<-read.csv("data/data_2019/data_used/Csv/EnteroD68_VP1.CSV")
HepatitisCvirus_1A<- read.csv("data/data_2019/data_used/Csv/HCV1A_CDS_Alignment.csv")
HepatitisCvirus_1B<- read.csv("data/data_2019/data_used/Csv/HCV1B_CDS_Alignment.csv")
HIV_1<- read.csv("data/data_2019/data_used/Csv/HIV1_FLT_2017_pol_DNA.csv")
#HumanBocavirus1_NS1<-read.csv("data/csv/HumanBocavirus1_NS1.CSV")
HumanBocavirus1_VP1<-read.csv("data/data_2019/data_used/Csv/HumanBoca_VP1.CSV")

humanparainfluenzavirus1_HN<-read.csv("data/data_2019/data_used/Csv/Parainfluenza1_HN.csv")
humanparainfluenzavirus3_HN<-read.csv("data/data_2019/data_used/Csv/Human_Parainfluenza3_HN.csv")
humanparainfluenzavirus1<-read.csv("data/data_2019/data_used/Csv/HumanParainfluenzavirus1_CDS_only.csv")
humanparainfluenzavirus3<-read.csv("data/data_2019/data_used/Csv/Parainfulenza3_CDS.csv")

InfluenzaAvirus_HA_H1N1<-read.csv("data/data_2019/data_used/Csv/InfluenzaAvirus_HA_H1N1_alignment.CSV")
InfluenzaAvirus_HA_H3N2<-read.csv("data/data_2019/data_used/Csv/InfluenzaAvirus_HA_H3N2_alignment.CSV")
InfluenzaAvirus_NA_H1N1<-read.csv("data/data_2019/data_used/Csv/InfluenzaAvirus_NA_H1N1_alignment.CSV")
InfluenzaAvirus_NA_H3N2<-read.csv("data/data_2019/data_used/Csv/InfluenzaAvirus_NA_H3N2_alignment.CSV")
InfluenzaBvirus_NA<-read.csv("data/data_2019/data_used/Csv/InfluenzaB_NA.CSV")
InfluenzaBvirus_HA<-read.csv("data/data_2019/data_used/Csv/InfluenzaB_HA.CSV")

HepatitisB_polymerase<-read.csv('data/data_2019/data_used/Csv/HBV_Polymerase.csv')
HepatitisB_precore<-read.csv("data/data_2019/data_used/Csv/HepatitisB_precore.csv")
HepatitisB_polymerase_truncated_precore<- read.csv("data/NESCENT_data/data_used/Csv/HepatitisB_polymerase_truncated_precore.csv")
HepatitisB_s<-read.csv('data/data_2019/data_used/Csv/HepatitisB_s.csv')
HepatitisB_pre_S<-read.csv('data/data_2019/data_used/Csv/HepatitisB_preS.csv')
HepatitisB_core<-read.csv("data/data_2019/data_used/Csv/HepatitisB_core.csv")

Humanherpesvirus2_glycoprotein_G<-read.csv("data/data_2019/data_used/Csv/Humanherpesvirus2_G.csv")
#Humanpapillomavirus16<-read.csv("data/csv/Humanpapillomavirus16.csv")
Humanpapillomavirus16_L1<-read.csv("data/data_2019/data_used/Csv/Humanpapilloma16_L1.csv")
Humanrespiratorysyncytialvirus<-read.csv("data/data_2019/data_used/Csv/HumanRespiratorySyncytial_CDS.csv") 
# Humanrespiratorysyncytialvirus_G<-read.csv("data/data_2019/data_used/Csv/Humanrespiratorysyncytialvirus_G.csv")
# JCpolyomavirus_VP1<-read.csv("data/csv/JCpolyomavirus_VP1.csv")
# Measles<-read.csv("data/csv/Measles.csv")
Measles_hemagglutinin_OR_haemagglutinin<-read.csv("data/data_2019/data_used/Csv/Measles_hemagglutinin_alignment.csv") 
ParvovirusB19_NS1<-read.csv("data/data_2019/data_used/Csv/ParovirusB19_NS1_alignment.csv")
ParvovirusB19_VP1<-read.csv("data/data_2019/data_used/Csv/Parovovirus_B19_VP1_alignment.csv") 
RhinovirusB<-read.csv("data/data_2019/data_used/Csv/RhinovirusB alignment.csv") 
#RhinovirusB_polyprotein<-read.csv("data/csv/RhinovirusB_polyprotein.csv")
RhinovirusC<-read.csv("data/data_2019/data_used/Csv/RhinovirusC.csv")
RotavirusA_VP6<-read.csv("data/data_2019/data_used/Csv/RotavirusA_VP6_alignment.csv")

Virus_info<- read.csv("data/list/Final_CpG_list.csv")
Virus_info$total=Virus_info$Number_of_Sequences*Virus_info$Number_of_Nucleotides

my.list <- list(DengueVirus1, DengueVirus2, DengueVirus3, DengueVirus4, HepatitisCvirus_1A, HepatitisCvirus_1B, HIV_1, humanparainfluenzavirus1_HN, humanparainfluenzavirus3_HN, humanparainfluenzavirus1, humanparainfluenzavirus3, InfluenzaAvirus_HA_H1N1,InfluenzaAvirus_HA_H3N2, InfluenzaAvirus_NA_H1N1, InfluenzaAvirus_NA_H3N2,InfluenzaBvirus_HA, InfluenzaBvirus_NA, EnterovirusA_VP1, EnterovirusA_VP2,EnterovirusB_VP1, EnterovirusB_VP2,EnterovirusC_VP1,EnterovirusC_VP2,EnterovirusD_VP1, Humanrespiratorysyncytialvirus, Measles_hemagglutinin_OR_haemagglutinin, RhinovirusB, RhinovirusC, RotavirusA_VP6, BKpolyomavirus_VP1, HumanBocavirus1_VP1, HepatitisB_polymerase,HepatitisB_precore,HepatitisB_polymerase_truncated_precore,HepatitisB_s,HepatitisB_pre_S,HepatitisB_core, Humanherpesvirus2_glycoprotein_G, Humanpapillomavirus16_L1, ParvovirusB19_NS1, ParvovirusB19_VP1)
#name.list <- list('DengueVirus1', 'DengueVirus2', 'DengueVirus3', 'DengueVirus4', 'humanparainfluenzavirus1_F', 'humanparainfluenzavirus1_HN', 'humanparainfluenzavirus3_HN', 'humanparainfluenzavirus1', 'humanparainfluenzavirus3', 'InfluenzaAvirus_HA_H1N1','InfluenzaAvirus_HA_H3N2', 'InfluenzaAvirus_NA_H1N1', 'InfluenzaAvirus_NA_H3N2','InfluenzaBvirus_HA', 'InfluenzaBvirus_NA', 'EnterovirusA_VP1', 'EnterovirusA_VP2','EnterovirusB_VP1', 'EnterovirusB_VP2','EnterovirusC_VP1','EnterovirusC_VP2','EnterovirusD_VP1', 'BKpolyomavirus_VP1', 'HumanBocavirus1_VP1', 'HepatitisB_polymerase','HepatitisB_precore','HepatitisB_ptp', 'HepatitisB_s','HepatitisB_pre_S','HepatitisB_core', 'Humanherpesvirus2_glycoprotein_G', 'Humanpapillomavirus16_L1', 'Humanrespiratorysyncytialvirus', 'Humanrespiratorysyncytialvirus_G', 'JCpolyomavirus_VP1', 'Measles_hemagglutin', 'ParvovirusB19_NS1', 'ParvovirusB19_VP1', 'RhinovirusB', 'RhinovirusB_polyprotein', 'RhinovirusC', 'RotavirusA_VP6')
name.list <- list('Dengue 1', 'Dengue 2', 'Dengue 3', 'Dengue 4', "HCV 1 A", "HCV 1 B", "HIV pol Gene", 'Human Parainfluenza 1 HN', 'Human Parainfluenza 3 HN', 'Human Parainfluenza 1', 'Human Parainfluenza 3','Influenza A HA H1N1','Influenza A HA H3N2', 'Influenza A NA H1N1', 'Influenza A NA H3N2','Influenza B HA', 'Influenza B NA', 'Entero A VP1 EVA71', 'Entero A VP2 EVA71','Entero B VP1 Echovirus30', 'Entero B VP2 Echovirus30','Entero C VP1 Polio2','Entero C VP2 Polio2','Entero D68 VP1', 'Human Respiratory Syncytial', 'Measles HH', 'Rhino B', 'Rhino C', 'Rota A VP6', 'Bk Polyoma VP1', 'Human Boca 1 VP1', 'Hepatitis B Polymerase','Hepatitis B Precore','Hepatitis B PTP','Hepatitis B S','Hepatitis B PreS','Hepatitis B Core', 'Herpes 2 Glycoprotein G', 'Human Papilloma 16 L1', 'Parvo B19 NS1', 'Parvo B19 VP1')
#Virus_info$total[33]=Virus_info$total[48] # influenza short
#Virus_info$total[36]=Virus_info$total[49]
reorg<-Virus_info[match(name.list, Virus_info$name),]
#reorg[!reorg$used == "no", ]
#reorg[1:9] <- NULL
#reorg[2:8] <- NULL
reorg$count=c(1:nrow(reorg))
data_points = data.frame("Count"= 1:length(my.list), "Virus"= 1:length(my.list))
count = 1

mean_or_median<- "mean"
for (data in my.list){  
  cpg.y<-subset(data, makesCpG==1)
  cpg.n<-subset(data, makesCpG==0)
  #subset further into letters nuclotideCpgforming or nucotideNonGpg
  AC<-subset(cpg.y, wtnt_consensus=='a')
  ANC<-subset(cpg.n, wtnt_consensus=='a')
  TC<-subset(cpg.y, wtnt_consensus=='t')
  TNC<-subset(cpg.n, wtnt_consensus=='t')
  
  
  #Function to help create errorbars
  sem<-function(x){
    return(sd(x,na.rm = FALSE)/sqrt(length(x)))
  }
  if(mean_or_median == "mean"){
    stats= mean
    error_bar= sem
  }
  if(mean_or_median== "median"){
    stats= median
    error_bar= sd
  }
  
  #making the data frames with all information about a, t, c, g 
  AllA = rbind(AC, ANC)
  AllT = rbind(TC, TNC)
  
  # for loops to caculate mean ans errorbars
  for (i in 1:length(AllA$makesCpG)) {
    if (AllA$makesCpG[i] == 1 && AllA$TypeOfSite[i] == "syn") {
      AllA_mean_value_syn_CpG <- stats(AllA$Freq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") )])
      AllA_sem_vals_syn_CpG<-error_bar(AllA$Freq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") )])

          }
    if (AllA$makesCpG[i] == 1 && AllA$TypeOfSite[i] == "nonsyn") {
      AllA_mean_value_nonsyn_CpG<- stats(AllA$Freq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "nonsyn") )])
      AllA_sem_vals_nonsyn_CpG<-error_bar(AllA$Freq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "nonsyn") )])
    }
    if (AllA$makesCpG[i] == 0 && AllA$TypeOfSite[i] == "syn") {
      AllA_mean_value_syn_nCpG <- stats(AllA$Freq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") )])
      
      AllA_sem_vals_syn_nCpG<-error_bar(AllA$Freq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") )])
    }
    if (AllA$makesCpG[i] == 0 && AllA$TypeOfSite[i] == "nonsyn") {
      AllA_mean_value_nonsyn_nCpG<- stats(AllA$Freq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "nonsyn") )])
      AllA_sem_vals_nonsyn_nCpG<-error_bar(AllA$Freq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "nonsyn") )])
    }
  }
  
  for (i in 1:length(AllT$makesCpG)) {
    if (AllT$makesCpG[i] == 1 && AllT$TypeOfSite[i] == "syn") {
      AllT_mean_value_syn_CpG <- stats(AllT$Freq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )])
      AllT_sem_vals_syn_CpG<-error_bar(AllT$Freq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )])
    }
    if (AllT$makesCpG[i] == 1 && AllT$TypeOfSite[i] == "nonsyn") {
      AllT_mean_value_nonsyn_CpG <- stats(AllT$Freq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "nonsyn") )])
      AllT_sem_vals_nonsyn_CpG<-error_bar(AllT$Freq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "nonsyn") )])
    }
    if (AllT$makesCpG[i] == 0 && AllT$TypeOfSite[i] == "syn") {
      AllT_mean_value_syn_nCpG<- stats(AllT$Freq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )])
      AllT_sem_vals_syn_nCpG<-error_bar(AllT$Freq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )])
    }
    if (AllT$makesCpG[i] == 0 && AllT$TypeOfSite[i] == "nonsyn") {
      AllT_mean_value_nonsyn_nCpG <- stats(AllT$Freq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "nonsyn") )])
      AllT_sem_vals_nonsyn_nCpG<-error_bar(AllT$Freq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "nonsyn") )])
    }
  }
  
  # There are the upper and lower limits of the error bar
  # AllA_LCLS = AllA$mean_value - AllA$sem_vals
  # AllA_UCLS = AllA$mean_value + AllA$sem_vals
  # 
  # AllT_LCLS = AllT$mean_value - AllT$sem_vals
  # AllT_UCLS = AllT$mean_value + AllT$sem_vals
  # 
  
  data_points$Virus[count]= name.list[count]
  data_points$AsynNC_C[count]= AllA_mean_value_syn_nCpG/AllA_mean_value_syn_CpG
  data_points$AnonsynNC_C[count]= AllA_mean_value_nonsyn_nCpG/AllA_mean_value_nonsyn_CpG
  
  data_points$TsynNC_C[count] = AllT_mean_value_syn_nCpG/AllT_mean_value_syn_CpG
  data_points$TnonsynNC_C[count] = AllT_mean_value_nonsyn_nCpG/AllT_mean_value_nonsyn_CpG
  
  
  data_points$AsynC_LCLS[count]= AllA_mean_value_syn_CpG - AllA_sem_vals_syn_CpG
  data_points$AnonsynC_LCLS[count]= AllA_mean_value_nonsyn_CpG - AllA_sem_vals_nonsyn_CpG
  data_points$AsynNC_LCLS[count]= AllA_mean_value_syn_nCpG - AllA_sem_vals_syn_nCpG
  data_points$AnonsynNC_LCLS[count]= AllA_mean_value_nonsyn_nCpG - AllA_sem_vals_nonsyn_nCpG
  data_points$TsynC_LCLS[count] = AllT_mean_value_syn_CpG - AllT_sem_vals_syn_CpG
  data_points$TnonsynC_LCLS[count] =AllT_mean_value_nonsyn_CpG - AllT_sem_vals_nonsyn_CpG
  data_points$TsynNC_LCLS[count] = AllT_mean_value_syn_nCpG - AllT_sem_vals_syn_nCpG 
  data_points$TnonsynNC_LCLS[count] =AllT_mean_value_nonsyn_nCpG - AllT_sem_vals_nonsyn_nCpG
  
  data_points$AsynC_UCLS[count]= AllA_mean_value_syn_CpG + AllA_sem_vals_syn_CpG
  data_points$AnonsynC_UCLS[count]= AllA_mean_value_nonsyn_CpG + AllA_sem_vals_nonsyn_CpG
  data_points$AsynNC_UCLS[count]= AllA_mean_value_syn_nCpG + AllA_sem_vals_syn_nCpG
  data_points$AnonsynNC_UCLS[count]= AllA_mean_value_nonsyn_nCpG + AllA_sem_vals_nonsyn_nCpG
  data_points$TsynC_UCLS[count] = AllT_mean_value_syn_CpG + AllT_sem_vals_syn_CpG
  data_points$TnonsynC_UCLS[count] =AllT_mean_value_nonsyn_CpG + AllT_sem_vals_nonsyn_CpG
  data_points$TsynNC_UCLS[count] = AllT_mean_value_syn_nCpG + AllT_sem_vals_syn_nCpG 
  data_points$TnonsynNC_UCLS[count] =AllT_mean_value_nonsyn_nCpG + AllT_sem_vals_nonsyn_nCpG
  count = count +1
}

#making new column of average fitness costs
data_points$averagefit<-0
for(i in 1:nrow(data_points)){
  avg <-c(data_points$AsynNC_C[i], data_points$AnonsynNC_C[i], data_points$TsynNC_C[i], data_points$TnonsynNC_C[i])
  for(inf in avg){
  }
  if(data_points$AsynNC_C[i] != "Inf" & data_points$AsynNC_C[i] !='NaN'){
    avgi<- c(data_points$AsynNC_C[i])
    #avg <-c(data_points$AnonsynNC_C[i], data_points$TsynNC_C[i], data_points$TnonsynNC_C[i])
  }
  if(data_points$AnonsynNC_C[i] != "Inf" & data_points$AnonsynNC_C[i] !='NaN'){
    avgi<-append(avgi, data_points$AnonsynNC_C[i])
  }
  if(data_points$TsynNC_C[i] != "Inf" & data_points$TsynNC_C[i] !='NaN'){
    avgi<- append(avgi,data_points$TsynNC_C[i])
  }
  if(data_points$TnonsynNC_C[i] != "Inf" & data_points$TnonsynNC_C[i] !='NaN'){
    avgi<- append(avgi,data_points$TnonsynNC_C[i])
  }
  
  data_points$averagefit[i]= mean(avgi)}

#sort average fitness costs to graph viruses in order of decreasing fitness
#data_points <-(data_points[order(data_points$averagefit,decreasing = TRUE),])
#rewrite count column
#data_points$Count <- 1:nrow(data_points)

#one of the error bars is too large at 2.47 e13 so we are placing it lower, but noting it high amount  
#data_points[20,12] = data_points[20,14]/300
#order by protein

##### protein test not used
# data_points$protein<-"missing"
# data_points$multNS = 0
# data_points$nucl = 0
# data_points$seq = 0
# #Virus_info<- read.csv("data/CpG_List.csv",na.strings = "NNN")
# for (i in 1:nrow(data_points)) {
#   treat=unlist(data_points[i,2])
#   
#   for (j in 1:nrow(Virus_info)){
#   # splitnameAll<-unlist(strsplit(as.character(Virus_info$name[j]),".fasta"))
#   # splitname<-unlist(strsplit(as.character(splitnameAll[1]),"_"))
#   truename<-Virus_info$nice_name[j]
#   #print(truename)
#   
#   if(treat==truename){
#     data_points$protein[i]<-as.character(Virus_info$gene[j])
#     data_points$nucl[i] = Virus_info$NucleotideNumber[j]
#     data_points$seq[i] = Virus_info$SeqNumber[j]
#     data_points$multNS[i]<-reorg$total[i]
#     data_points$multNS[i]<-reorg$total[i]
#    
#   }
#   } 
# }
df <- apply(data_points,2,as.character)
write.csv(df, file = "output/All_Data/Costly/alldatapoints.csv")

#data_points$Virus<-c('Dengue1', 'Dengue2', 'Dengue3', 'Dengue4', 'humanparainfluenza1_F', 'humanparainfluenza1_HN', 'humanparainfluenza3_HN', 'InfluenzaA_HA_H1N1','InfluenzaA_HA_H3N2', 'InfluenzaA_NA_H1N1', 'InfluenzaA_NA_H3N2','InfluenzaB_HA', 'InfluenzaB_NA', 'EnteroA_VP1', 'EnteroA_VP2','EnteroB_VP1', 'EnteroB_VP2','EnteroC_VP1','EnteroC_VP2','EnteroD_VP1', 'BK_polyoma_VP1', 'HumanBoca1_VP1', 'HepatitisB_polymerase','HepatitisB_precore','HepatitisB_polymerase_truncated_precore','HepatitisB_s','HepatitisB_pre_S','HepatitisB_core', 'Humanherpes2_glycoprotein_G', 'Humanpapilloma16_L1', 'Humanrespiratorysyncytial', 'Humanrespiratorysyncytial_G', 'JCpolyoma_VP1', 'Measles_hemagglutin', 'ParvoB19_NS1', 'ParvoB19_VP1', 'RhinoB', 'RhinoB_polyprotein', 'RhinoC', 'RotaA_VP6', 'humanparainfluenza1', 'humanparainfluenza3')

#data_points <-(data_points[order(data_points$protein),])
#rewrite count column
#data_points$Count <- 1:nrow(data_points)


#print(data_points$TnonsynNC_LCLS/data_points$TnonsynC_LCLS)
# graphing 
png("output/All_data/Costly/Costly_Graph_AllR_9_25.png", width = 15, height = 8, units = "in", res= 500)
#--------------------
par(mar=c(0,2,3,2), oma=c(6,4,1,1), mfrow=c(2,1))#, bg = "darkseagreen1"
#changed mar(0,2,3,2) oma(6,4,1,1)
x <- data_points$Count# test data 
y <-  data_points$AsynNC_C
plot(x,y, type = "n", log ='y' ,main="Cost of CpG-Creating Mutations", xlab=" ", yaxt = "n", ylab="Costly", xaxt = "n", ylim=c(0.3, 200), xlim=c(2, length(my.list) +5.5), las= 1, cex.main=3) 

abline(v=c(1.5,2.5,3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,12.5,13.5,14.5,16.5,17.5,18.5,19.5,15.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5), col="grey", lty=c(1))


u <- par('ylog') 
rect(-1.05, .0001, 4.5, 1570, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(6.5, .0001, 7.5, 1570, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(11.5, .0001, 17.5, 1570, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(24.5, .0001, 25.5, 1570, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(26.5, .0001, 28.5, 1570, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(29.5, .0001, 30.5, 1570, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(31.5, .0001, 37.5, 1570, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(38.5, .0001, 39.5, 1570, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(41.5, .0001, 50.5, 1570, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))

# aty <- axTicks(2)
# aty <- axTicks(2)
# labels <- sapply(aty,function(i)
#   as.expression(bquote(10^ .(i)))
# )
# axis(2,at=aty,labels=labels)
points(data_points$Count-.3, data_points$AsynNC_C, col='red', pch= 19)
#points(data_points$Count -.1, data_points$AnonsynNC_C, col= "green", pch=19)
points(data_points$Count+.1, data_points$TsynNC_C, col= "blue", pch=19)
#points(data_points$Count + .3, data_points$TnonsynNC_C, col= "purple", pch=19)
# par(new=TRUE)

# hack: we draw arrows but with very special "arrowheads"
arrows(data_points$Count-.3, data_points$AsynNC_LCLS/data_points$AsynC_LCLS, data_points$Count -.3, data_points$AsynNC_UCLS/data_points$AsynC_UCLS, length=0.05, angle=90, code=3, col= "red")
#arrows(data_points$Count-.1, data_points$AnonsynNC_LCLS/data_points$AnonsynC_LCLS, data_points$Count -.1, data_points$AnonsynNC_UCLS/data_points$AnonsynC_UCLS, length=0.05, angle=90, code=3, col= "green")
arrows(data_points$Count+.1, data_points$TsynNC_LCLS/data_points$TsynC_LCLS, data_points$Count+.1, data_points$TsynNC_UCLS/data_points$TsynC_UCLS, length=0.05, angle=90, code=3, col= "blue")
#arrows(data_points$Count+.3, data_points$TnonsynNC_LCLS/data_points$TnonsynC_LCLS, data_points$Count +.3, data_points$TnonsynNC_UCLS/data_points$TnonsynC_UCLS, length=0.05, angle=90, code=3, col= "purple")

axis(2, at = c(.01, 0.5,1,2,5,10,20,50,100,200), labels = c(0.01, 0.5,1,2,5,10,20,50,100,200),  las=2)
#axis.break(2, 0.007,breakcol="black",style="slash")
#mtext('No nonCpG \n  or  CpG     \n    mutations  ', side=2, line=.005, at=0.002, las=1.1, cex = .7)
#mtext('No nonCpG \n   mutations   ', side=2, line=.005, at=0.0005, las=1.1, cex = .7)

#axis.break(2, 200,breakcol="black",style="slash")
#mtext('2.47e+13', side=2, line=.005, at=300, las=1.1, cex = .9)

mtext('Relative Cost of CpG \n Creating Mutations', side=2, line=3, at=5, las=0, cex = 1.3)

#mtext('No CpG \n mutations ', side=2, line=.005, at=700, las=1.1, cex = .7)

abline(h=c(0.002, 0.0005, .01,0.5,1,2,5,10,20,50,100, 200), col="grey", lty=c(2,2))

abline(v=c(1.5,2.5,3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,12.5,13.5,14.5,16.5,17.5,18.5,19.5,15.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5), col="grey", lty=c(1))
#abline(v=c(6.5,10.5,14.5,20.5, 27.5, 41.5), col = "darkgreen", lwd=2)
abline(v=c(4.5,6.5,7.5,11.5,17.5,24.5,25.5,26.5,28.5,29.5,30.5,31.5,37.5,38.5,39.5,41.5), col = "bisque4", lwd=2)


abline(h = 1, col ="darkslategrey", lwd = 2)
#rect(0,800, 800,800, col = "darkseagreen1")
# xlab="Virus "
text(45, 2, " CpG Mutation \n More Costly", cex = 1, font = 2)
text(42.5, 2, "↑", cex = 3, font = 2)
text(45, .5, " CpG Mutation \n Less Costly", cex = 1, font = 2)
text(42.5, .5, "↓", cex = 3, font = 2)
#axis(1, at=1:length(my.list), labels=data_points$Virus, las= 2)
legend('topright', legend=c("A -> G Syn","", "T -> C Syn", " "),
       col=c("red", "white", "blue", "white"), lty=1, lwd= 3, cex = 1, pt.cex = 999)


par(mar=c(10,2,0,2))
plot(reorg$count,reorg$total, ylim=c(60000, 45000000), type = "n", log ='y' , xlab=" ", yaxt = "n",  xaxt = "n", xlim=c(2, length(my.list) +5.5), las= 1, cex.main=3) 
mtext('Amount of \n Data', side=2, line=3, at=1100000, las=0, cex = 1.3)

abline(v=c(1.5,2.5,3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,12.5,13.5,14.5,16.5,17.5,18.5,19.5,15.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5), col="grey", lty=c(1))
abline(v=c(4.5,6.5,7.5,11.5,17.5,24.5,25.5,26.5,28.5,29.5,30.5,31.5,37.5,38.5,39.5), col = "bisque4", lwd=2)

abline(h = 55988780, col ="darkslategrey", lwd = 2, lty = 2)
axis(2, at = c(60000,150000,1000000,5000000,25000000), labels = c("60,000","150k","1 mil","5 mil","25 mil"),  las=2)
abline(h=c(60000,150000,1000000,5000000,25000000), col="grey", lty=c(2,2))

u <- par('ylog') 

rect(-1.05, .0001, 4.5, 558988780, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(6.5, .0001, 7.5, 558988780, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(11.5, .0001, 17.5, 558988780, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(24.5, .0001, 25.5, 558988780, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(26.5, .0001, 28.5, 558988780, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(29.5, .0001, 30.5, 558988780, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(31.5, .0001, 37.5, 558988780, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(38.5, .0001, 39.5, 558988780, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(41.5, .0001, 50.5, 558988780, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
points(reorg$count,reorg$total, pch = 16)

axis(1, at=1:length(my.list), labels=data_points$Virus, las= 2)

dev.off()
# add horsontal gray lines

