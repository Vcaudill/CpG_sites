
#This function is going to read the data from the csv files 
#makesCpG <-> ancestor_makesCpG
#wtnt_consensus<->ancestor
#TypeOfSite<->ancestor_TypeOfSite
#Freq <-> qFreq
Tables = function(truename, csv_file){ 
  #setwd("..")
  #truename ="Humanherpesvirus2_gD"
  truenamecsv= paste(truename, ".csv", sep="")
  print(truenamecsv)
  DF<- read.csv(csv_file)
  #if(DF<- read.csv(paste(place, truenamecsv, sep="")) == FALSE) {
    #DF<- read.csv(paste("data/data_2019/Csv/new_for_costly", truenamecsv, sep=""))}
  #}
  #load data as DF
  # Wilcox test
  #Wilcox_test(DF, truename)
  return(DF)
  }


Wilcox_test = function(data, truename){
  
  #set output pdf file name
  library(graphics)
  library(plyr)
  library(dplyr)
 
  pVals = c()
  shrtval = 0
  options(scipen=999)
  #prevents pvalues from becoming scientific notation. 
  
  array1 = data$aFreq[data$ancestor =="a" & data$ancestor_TypeOfSite == 'syn' & data$ancestor_makesCpG == 1]
  array2 = data$aFreq[data$ancestor =="a" & data$ancestor_TypeOfSite == 'syn' & data$ancestor_makesCpG == 0]
  array3 = data$aFreq[data$ancestor =="a" & data$ancestor_TypeOfSite == 'nonsyn' & data$ancestor_makesCpG == 1]
  array4 = data$aFreq[data$ancestor =="a" & data$ancestor_TypeOfSite == 'nonsyn' & data$ancestor_makesCpG == 0]
  syna = data$aFreq[data$ancestor =="a" & data$ancestor_TypeOfSite == 'syn']
  nonsyna = data$aFreq[data$ancestor =="a" & data$ancestor_TypeOfSite == 'nonsyn']
  CpGa = data$aFreq[data$ancestor =="a" & data$ancestor_makesCpG == 1]
  nonCpGa = data$aFreq[data$ancestor =="a" &  data$ancestor_makesCpG == 0]
  
  print("For a: Comparing makes CpG with noCpG (syn). Wilcox test less: red/blue")
  print(wilcox.test(array1, array2, alternative='less'))
  pVals = c(pVals,format(wilcox.test(array1, array2, alternative='less')$p.value, nsmall = 6))
  print(pVals)
  print("For a: Comparing CpG with noCpG (nonsyn). Wilcox test less: yellow/green")
  print(wilcox.test(array3, array4, alternative='less'))
  pVals = c(pVals,format(wilcox.test(array3, array4, alternative='less')$p.value, nsmall = 6))
  print(pVals)
  print("For a: Comparing  syn to nonsyn. Wilcox test greater red&blue vs yellow&green")
  print(wilcox.test(syna, nonsyna, alternative='greater'))
  print(wilcox.test(syna, nonsyna, alternative='greater')$p.value)
  pVals = c(pVals,format(wilcox.test(syna, nonsyna, alternative='greater')$p.value, nsmall = 6))
  print(pVals)
  
  array5 = data$aFreq[data$ancestor =="t" & data$ancestor_TypeOfSite == 'syn' & data$ancestor_makesCpG == 1]
  array6 = data$aFreq[data$ancestor =="t" & data$ancestor_TypeOfSite == 'syn' & data$ancestor_makesCpG == 0]
  array7 = data$aFreq[data$ancestor =="t" & data$ancestor_TypeOfSite == 'nonsyn' & data$ancestor_makesCpG == 1]
  array8 = data$aFreq[data$ancestor =="t" & data$ancestor_TypeOfSite == 'nonsyn' & data$ancestor_makesCpG == 0]
  synt = data$aFreq[data$ancestor =="t" & data$ancestor_TypeOfSite == 'syn']
  nonsynt = data$aFreq[data$ancestor =="t" & data$ancestor_TypeOfSite == 'nonsyn']
  CpGt = data$aFreq[data$ancestor =="t" & data$ancestor_makesCpG == 1]
  nonCpGt = data$aFreq[data$ancestor =="t" &  data$ancestor_makesCpG == 0]
  
  print("For t: Comparing makes CpG with noCpG (syn). Wilcox test less: red/blue")
  print(wilcox.test(array5, array6, alternative='less'))
  pVals = c(pVals,format(wilcox.test(array5, array6, alternative='less')$p.value, nsmall = 6))
  print("For t: Comparing CpG with noCpG (nonsyn). Wilcox test less: yellow/green")
  print(wilcox.test(array7, array8, alternative='less'))
  pVals = c(pVals,format(wilcox.test(array7, array8, alternative='less')$p.value, nsmall = 6))
  print("For t: Comparing  syn to nonsyn. Wilcox test greater red&blue vs yellow&green")
  print(wilcox.test(synt, nonsynt, alternative='greater'))
  pVals = c(pVals,format(wilcox.test(synt, nonsynt, alternative='greater')$p.value, nsmall = 6))
  
  Pvalues= c(pVals)
  #save Pvalues into list
  
  return(Pvalues)
  #makeTable(Pvalues, truename)
  
}

makeTable <- function(Pvalues, truename, nice_name,output){
  options(scipen = 999)
  #setwd("output/redeploy/")
  truenamepdf= paste(output,truename,"_ancestor.pdf",sep="")
  truenamepng= paste(truename,"_ancestor_tables", ".png", sep="")
  #print(truenamepdf)
  #prevents pvalues from becoming scientific notation
  options(warn=-1)
  #suppress warnings
  
  #setwd("~/Desktop/Something_Cool-CpG_Sites-/Tables")
  #table construct
  pdf(truenamepdf, width = 7, height= 5)
  #png(truenamepng, width = 6.75, height = 6.75, units = "in", res= 300)
  col1 <- c("A->G", "T->C")
  col2 <- c("Syn: CpG v NonCpG", "NonSyn: CpG v NonCpG", "Syn v NonSyn")
  ycoor <- c(4*100/5+.7 , 3*100/5 + 5.1, 3*100/5 - 10, 2*100/5 -6.4, 1*100/5-1.3, 100/5-14- 2.9, 100)
  ycoorb <- c(4*100/5+.7 , 3*100/5 + 5.1, 3*100/5 - 10.6, 2*100/5 -6.4, 1*100/5-1.3, 100/5-14- 2.9, 100)
  df = data.frame(col1, col2, Pvalues)
  
  par(xpd=F)
  plot(1, 2, xlim=c(0,100),ylim=c(0,100), col=0, xaxt="n", yaxt="n", xlab="", ylab="")
  title(main = nice_name, family = "Times", adj = 0.5, cex.main= 2)
  abline(v = 100/5)
  abline(v = 2*100/3)
  abline(h = 100-100/7 + 3)
  abline(h = 100- 2*100/7+2)
  abline(h = 100 - 3*100/7+1)
  abline(h= 100 - 4*100/7 -1)
  abline(h= 100 -5*100/7 - 2)
  abline(h= 100 -6*100/7 - 3)
  
  
  text(x=100/7- 6, y= 5*100/5-3, "Mutation Type")
  text(x=3*100/7, y = 5*100/5-3, "Comparison")
  text(x= 6*100/7, y=5*100/5-3, "P-Value")
  rect(xleft = -4, xright = 100/5, ybottom =42, ytop =100-100/7+3 , col = "white")
  text(x= 100/12, y= 3*100/5+5, "A->G", cex = 1.7, family = "Times")
  rect(xleft = -4, xright = 100/5, ybottom =-4, ytop =42 , col = "white")
  text(x= 100/12, y = 1*100/5 - 1, "T->C", cex = 1.7, family ='Times')
  text(x= 3*100/7, y = ycoor[1], labels= col2[1])
  text(x= 3*100/7, y = ycoor[2], labels= col2[2])
  text(x= 3*100/7, y = ycoor[3], labels= col2[3])
  text(x= 3*100/7, y = ycoor[4], labels= col2[1])
  text(x= 3*100/7, y = ycoor[5], labels= col2[2])
  text(x= 3*100/7, y = ycoor[6], labels= col2[3])
  
  num <- 1
  for (i in Pvalues){
    #i = format(i, nsmall = 6)
    
    
    library(scales)
    if (i < 0.01){
      a = 0.45
      i = "< 0.01"
    }
    else if(i <0.05){
      a = 0.25
      i = as.numeric(i)
      i = signif(i, digits = 3)
    }
    else if(i >0.05){
      a = 0.1
      i = as.numeric(i)
      i = signif(i, digits = 3)
      print(i)
    }
    
    rect(xleft = 2*100/3, xright = 200, ybottom = ycoorb[num]-7.3, ytop = ycoor[num]+8, col = alpha("deepskyblue1", a), border = col)
    text(x= 6*100/7, y =ycoor[num], labels = i)
    num = num + 1 
  }
  print("end")
  #dev.copy(pdf, truenamepng)
  dev.off()
  
}

#loop through namelist (all viruses)
# 
# hyphy_virus<-read.csv("data/list/old_lists/CpG_List_NewdataRW.csv")
# for(i in 1:nrow(hyphy_virus)){
#   nice_name <- as.character(hyphy_virus$nice_name[i])
#   print(hyphy_virus$name[i])
#   name = as.character(hyphy_virus$name[i])
#   #if(grep(".fasta",hyphy_virus$name[i])){
#   #  splitname<-unlist(strsplit(as.character(hyphy_virus$name[i]),".fasta"))
#   truename<- name
#    # }
#  # else{(truename<-name)}
#   print(truename)
#   if (truename == "Humanherpesvirus2_gD") {
#     next
#     }
# #     
#   place= hyphy_virus$place[i]
#   DF=Tables(truename, place)
#   Pvalues=Wilcox_test(DF, truename)
#   output = hyphy_virus$output_path[i]
#   makeTable(Pvalues, truename, nice_name, output)
#   data_place = hyphy_virus$place[i]
#   data_output = paste(output, "M_frequency_graphs/", sep="")
#   source(("R_scripts/graphs/M_frequency_graph.R"))
#   comparing_CpG_Syn_Nonsyn_new(truename,nice_name,data_place,data_output)
#   }
# readfile <- function(truenamecsv){
#   tryCatch(
#     # This is what I want to do...
#     {
#       DF<- read.csv(paste("new_data/Consensus/", truenamecsv, sep=""))
#       return(DF)
#     },
#     # ... but if an error occurs, tell me what happened: 
#     error=function(error_message) {
#       message("This is my custom message.")
#       message("Error in file(file, "rt") : cannot open the connection")
#       message(error_message)
#       return(NA)
#     }
#   )

