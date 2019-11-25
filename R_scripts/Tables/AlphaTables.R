#these are updated name lists
#my.list <- list(DengueVirus1, DengueVirus2, DengueVirus3, DengueVirus4, HepatitisCvirus_1A, HepatitisCvirus_1B, HIV_1, humanparainfluenzavirus1_F, humanparainfluenzavirus1_HN, humanparainfluenzavirus3_HN, humanparainfluenzavirus1, humanparainfluenzavirus3, InfluenzaAvirus_HA_H1N1,InfluenzaAvirus_HA_H3N2, InfluenzaAvirus_NA_H1N1, InfluenzaAvirus_NA_H3N2,InfluenzaBvirus_HA, InfluenzaBvirus_NA, EnterovirusA_VP1, EnterovirusA_VP2,EnterovirusB_VP1, EnterovirusB_VP2,EnterovirusC_VP1,EnterovirusC_VP2,EnterovirusD_VP1, Humanrespiratorysyncytialvirus, Humanrespiratorysyncytialvirus_G, Measles_hemagglutinin_OR_haemagglutinin, RhinovirusB, RhinovirusC, RotavirusA_VP6, BKpolyomavirus_VP1, HumanBocavirus1_VP1, HepatitisB_polymerase,HepatitisB_precore,HepatitisB_polymerase_truncated_precore,HepatitisB_s,HepatitisB_pre_S,HepatitisB_core, Humanherpesvirus2_glycoprotein_G, Humanpapillomavirus16_L1, ParvovirusB19_NS1, ParvovirusB19_VP1)
#name.list <- list('Dengue 1', 'Dengue 2', 'Dengue 3', 'Dengue 4', "HCV1A", "HCV1B", "HIV pol gene", 'Human Parainfluenza 1 F', 'Human Parainfluenza 1 HN', 'Human Parainfluenza 3 HN', 'Human Parainfluenza 1', 'Human Parainfluenza 3', 'Influenza A H1N1 HA','Influenza A H3N2 HA', 'Influenza A H1N1 NA', 'Influenza A H3N2 NA','Influenza B HA', 'Influenza B NA', 'Entero A VP1', 'Entero A VP2','Entero B VP1', 'Entero B VP2','Entero C VP1','Entero C VP2','Entero D VP1', 'Human Respiratory Syncytial', 'Human Respiratory Syncytial G', 'Measles HH', 'Rhino B', 'Rhino C', 'Rota A VP6', 'Bk Polyoma VP1', 'Human Boca 1 VP1', 'Hepatitis B Polymerase','Hepatitis B Precore','Hepatitis B PTP','Hepatitis B S','Hepatitis B PreS','Hepatitis B Core', 'Human Herpes 2 glycoprotein G', 'Human Papilloma 16 L1', 'Parvo B19 NS1', 'Parvo B19 VP1')

######test something
# 
#path<-paste("Hyphy/","BK_stuff","/", "fubar","/datamonkey-table.csv",sep='')
#datamonkey<-read.csv(path)
#cpg.yca<-subset(datamonkey, makesCpG==1 & potential_CpG == "yes" & last_nuc == "a" )
#makes CpG in genetic code and is 2 fold 
#cpg.ynca<-subset(datamonkey, makesCpG==0 & potential_CpG == "yes" & last_nuc == "a")
#does not make CpG in genetic code 
#cpg.yct<-subset(datamonkey, makesCpG==1 & potential_CpG == "yes" & last_nuc == "t" )
# #makes CpG in genetic code and is 2 fold 
#cpg.ynct<-subset(datamonkey, makesCpG==0 & potential_CpG == "yes" & last_nuc == "t")
# 


Wilcox_test = function(hyphy_virus_name, datamonkey, cpg.yca, cpg.yct, cpg.ynca, cpg.ynct){
  
  #set output pdf file name
  library(graphics)
  library(plyr)
  library(dplyr)
  
  pVals = c()
  shrtval = 0
  options(scipen=999)
  #prevents pvalues from becoming scientific notation. 
  
  array1 = cpg.yca$X.alpha.
  array2 = cpg.ynca$X.alpha.
  array3 = cpg.yct$X.alpha.
  array4 = cpg.ynct$X.alpha.
  #testing array1 vs array2; array3 vs array4
  
  print("For a: Comparing makes CpG with noCpG. Wilcox test less: red/blue")
  print(wilcox.test(array1, array2, alternative='less'))
  pVals = c(pVals,format(wilcox.test(array1, array2, alternative='less')$p.value, nsmall = 6))
  print(pVals)
  print("For t: Comparing CpG with noCpG. Wilcox test less: yellow/green")
  print(wilcox.test(array3, array4, alternative='less'))
  pVals = c(pVals,format(wilcox.test(array3, array4, alternative='less')$p.value, nsmall = 6))
  print(pVals)

  Pvalues= c(pVals)
  #save Pvalues into list

  options(scipen = 999)
  
  #setwd("output/redeploy/")
  truenamepdf= paste("output/Datamonkey/AlphaTables/",hyphy_virus_name,".pdf",sep="")
  truenamepng= paste(hyphy_virus_name,"tables", ".png", sep="")
  #print(truenamepdf)
  #prevents pvalues from becoming scientific notation
  options(warn=-1)
  #suppress warnings
  
  #setwd("~/Desktop/Something_Cool-CpG_Sites-/Tables")
  #table construct
  pdf(truenamepdf, width = 7, height= 5)
  #png(truenamepng, width = 6.75, height = 6.75, units = "in", res= 300)
  col1 <- c("A->G", "T->C")
  #ycoor <- c(4*100/5+.7 , 3*100/5 + 5.1, 3*100/5 - 10, 2*100/5 -6.4, 1*100/5-1.3, 100/5-14- 2.9, 50)#7 values for 7 rows in wilcox tables
  ycoor <- c(4*100/5+.7 , 3*100/5 + 5.1, 3*100/5 - 10)
  ycoorb <- c(4*100/5+.7 , 3*100/5 + 5.1, 3*100/5 - 10.6, 2*100/5 -6.4, 1*100/5-1.3, 100/5-14- 2.9, 50)
  df = data.frame(col1, Pvalues)
  
  par(xpd=F)
  plot(1, 2, xlim=c(0,100),ylim=c(60,100), col=0, xaxt="n", yaxt="n", xlab="", ylab="")
  title(main = nice_name, family = "Times", adj = 0.5, cex.main= 2)
  abline(v = 100/2)
  
  abline(h = 100-100/7 + 3)
  #abline(h = 100- 2*100/7+2)
  #abline(h = 100 - 3*100/7+1)
  abline(h= 100 - 100/7 -12.5)
  #abline(h= 100 -5*100/7 - 2)
  #abline(h= 100 -6*100/7 - 3)
  
  
  text(x=100/7- 6, y= 5*100/5-3, "Mutation Type")
  text(x= 6*100/7, y=5*100/5-3, "P-Value")
  #rect(xleft = -4, xright = 100/5, ybottom =42, ytop =100-100/7+3 , col = "white")
  text(x= 100/12, y= 3*100/4+5, "A->G", cex = 1.2, family = "Times")
  #rect(xleft = -4, xright = 100/5, ybottom =-4, ytop =42 , col = "white")
  text(x= 100/12, y = 3*100/4 -7, "T->C", cex = 1.2, family ='Times')
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
      i = signif(i, digits = 4)
    }
    else if(i >0.05){
      a = 0.1
      i = as.numeric(i)
      i = signif(i, digits = 4)
      print(i)
    }
    
    rect(xleft = 120, xright = -20, ybottom = ycoorb[num]-7.3, ytop = ycoor[num]+8, col = alpha("deepskyblue1", a), border = col)
    text(x= 6*100/7, y =ycoor[num], labels = i)
    num = num + 1 
  }
  print("end")
  # pdf(truenamepdf, width = 7, height= 5)
  # dev.copy(pdf, truenamepng)
  dev.off()
  
}

#loop through namelist (all viruses)
# for(i in 1:nrow(Virus_info)){
#   #print(Virus_info$name[i])
#   name = as.character(Virus_info$name[i])
#   splitname<-unlist(strsplit(as.character(Virus_info$name[i]),".fasta"))
#   truename<-splitname[1]
#   print(truename)
#   if (truename == "Humanherpesvirus2_gD") {
#     next
#     }
#     
#   DF=Tables(truename)
#   Pvalues=Wilcox_test(DF, truename)
#   makeTable(Pvalues, truename)
#   }

