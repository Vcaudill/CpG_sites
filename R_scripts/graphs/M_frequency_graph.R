
comparing_CpG_Syn_Nonsyn_new = function(truename,nice_name,data_place,data_output){
  # place is srting of text like where file should be /new_data/folder
  data<- read.csv(paste(data_place, truename,".csv", sep=""))
  #subset into two groups yes makes cpg and no cpg
  cpg.y<-subset(data, makesCpG==1)
  cpg.n<-subset(data, makesCpG==0)
  #subset further into letters nuclotideCpgforming or nucotideNonGpg
  AC<-subset(cpg.y, wtnt_consensus=='a')
  ANC<-subset(cpg.n, wtnt_consensus=='a') 
  GC<-subset(cpg.y, wtnt_consensus=='g') 
  GNC<-subset(cpg.n, wtnt_consensus=='g')
  TC<-subset(cpg.y, wtnt_consensus=='t')
  TNC<-subset(cpg.n, wtnt_consensus=='t')
  CC<-subset(cpg.y, wtnt_consensus=='c')
  CNC<-subset(cpg.n, wtnt_consensus=='c')
  
  #Function to help create errorbars
  sem<-function(x){
    return(sd(x,na.rm = FALSE)/sqrt(length(x)))
  }
  mean_or_median = "mean"
  
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
  AllC = CNC
  AllG = GNC
  
  # added new colmuns for mean an error bars
  AllA$mean_value <- .01
  AllA$sem_vals<- 0
  AllT$mean_value <- .01
  AllT$sem_vals<- 0
  AllC$mean_value <- .01
  AllC$sem_vals<- 0
  AllG$mean_value <- .01
  AllG$sem_vals<- 0
  
  AllA$median_value <- .01
  AllA$sd_vals<- 0
  AllT$median_value <- .01
  AllT$sd_vals<- 0
  AllC$median_value <- .01
  AllC$sd_vals<- 0
  AllG$median_value <- .01
  AllG$sd_vals<- 0

  
  AllAT = rbind(AllA, AllT)
  AllATCG = rbind(AllA, AllT, AllC, AllG)
  
  for (i in 1:length(AllA$makesCpG)) {
    if (AllA$makesCpG[i] == 1 && AllA$TypeOfSite[i] == "syn") {
      AllA$graphit[i] <- 2
      AllA$mean_value[i] <- stats(AllA$Freq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") )])
      AllA$sem_vals[i]<-error_bar(AllA$Freq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") )])
    }
    if (AllA$makesCpG[i] == 1 && AllA$TypeOfSite[i] == "nonsyn") {
      AllA$graphit[i] <- 4
      AllA$mean_value[i] <- stats(AllA$Freq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "nonsyn") )])
      AllA$sem_vals[i]<-error_bar(AllA$Freq[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "nonsyn") )])
    }
    if (AllA$makesCpG[i] == 0 && AllA$TypeOfSite[i] == "syn") {
      AllA$graphit[i] <- 1
      AllA$mean_value[i] <- stats(AllA$Freq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") )])
      AllA$sem_vals[i]<-error_bar(AllA$Freq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") )])
    }
    if (AllA$makesCpG[i] == 0 && AllA$TypeOfSite[i] == "nonsyn") {
      AllA$graphit[i] <- 3
      AllA$mean_value[i] <- stats(AllA$Freq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "nonsyn") )])
      AllA$sem_vals[i]<-error_bar(AllA$Freq[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "nonsyn") )])
    }
  }
  
  for (i in 1:length(AllT$makesCpG)) {
    if (AllT$makesCpG[i] == 1 && AllT$TypeOfSite[i] == "syn") {
      AllT$graphit[i] <- 2
      AllT$mean_value[i] <- stats(AllT$Freq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )])
      AllT$sem_vals[i]<-error_bar(AllT$Freq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )])
    }
    if (AllT$makesCpG[i] == 1 && AllT$TypeOfSite[i] == "nonsyn") {
      AllT$graphit[i] <- 4
      AllT$mean_value[i] <- stats(AllT$Freq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "nonsyn") )])
      AllT$sem_vals[i]<-error_bar(AllT$Freq[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "nonsyn") )])
    }
    if (AllT$makesCpG[i] == 0 && AllT$TypeOfSite[i] == "syn") {
      AllT$graphit[i] <- 1
      AllT$mean_value[i] <- stats(AllT$Freq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )])
      AllT$sem_vals[i]<-error_bar(AllT$Freq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )])
    }
    if (AllT$makesCpG[i] == 0 && AllT$TypeOfSite[i] == "nonsyn") {
      AllT$graphit[i] <- 3
      AllT$mean_value[i] <- stats(AllT$Freq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "nonsyn") )])
      AllT$sem_vals[i]<-error_bar(AllT$Freq[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "nonsyn") )])
    }
  }
  
  for (i in 1:length(AllC$makesCpG)) {
    
    if (AllC$makesCpG[i] == 0 && AllC$TypeOfSite[i] == "syn") {
      AllC$graphit[i] <- 1
      AllC$mean_value[i] <- stats(AllC$Freq[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "syn") )])
      AllC$sem_vals[i]<-error_bar(AllC$Freq[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "syn") )])
    }
    if (AllC$makesCpG[i] == 0 && AllC$TypeOfSite[i] != "syn") {
      AllC$graphit[i] <- 3
      AllC$mean_value[i] <- stats(AllC$Freq[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "nonsyn") )])
      AllC$sem_vals[i]<-error_bar(AllC$Freq[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "nonsyn") )])
    }
  }
  
  for (i in 1:length(AllG$makesCpG)) {
    
    if (AllG$makesCpG[i] == 0 && AllG$TypeOfSite[i] == "syn") {
      AllG$graphit[i] <- 1
      AllG$mean_value[i] <- stats(AllG$Freq[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "syn") )])
      AllG$sem_vals[i]<-error_bar(AllG$Freq[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "syn") )])
    }
    if (AllG$makesCpG[i] == 0 && AllG$TypeOfSite[i] != "syn") {
      AllG$graphit[i] <- 3
      AllG$mean_value[i] <- stats(AllG$Freq[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "nonsyn") )])
      AllG$sem_vals[i]<-error_bar(AllG$Freq[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "nonsyn") )])
    }
  }
  
  AllA$LCLS = AllA$mean_value - AllA$sem_vals
  AllA$UCLS = AllA$mean_value + AllA$sem_vals
  
  AllT$LCLS = AllT$mean_value - AllT$sem_vals
  AllT$UCLS = AllT$mean_value + AllT$sem_vals
  
  AllC$LCLS = AllC$mean_value - AllC$sem_vals
  AllC$UCLS = AllC$mean_value + AllC$sem_vals
  
  AllG$LCLS = AllG$mean_value - AllG$sem_vals
  AllG$UCLS = AllG$mean_value + AllG$sem_vals
  ####################################################################################
  #layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
  # make 2 dataset 1 conting everything except 0's (Freq needs to be log)
  # datea set just for Freq = 0's
  # graph it with breaks 
  # add a key
  # make it more like jitter
  # try to make points transparent 30% via alpha or make error darker in color
  # make bars the same width
  # try to get good png size 
  
  library(scales)
  library(plotrix)
  library(sfsmisc)
  
  truenamepng = paste(data_output,truename,".png",sep="")
  png(truenamepng, width = 6.75, height = 6.75, units = "in", res= 300)
  par(mfrow=c(2,2), mar=c(4.1, 4.1, 1.9, 0.8),oma=c(0.1,0.1,1.5,0.1)) 
  palette(alpha(c("#99FF99","#9999FF","#FF9900","#FF3300"),0.3))

  
  plot(jitter(AllA$graphit),AllA$Freq + 0.0001,log='y',col=factor(AllA$graphit),pch=16, main="A->G",xlab = " ", ylab = "Mutation Frequency", yaxt="n", xaxt="n", ylim=c(0.0001, 0.5))

  points(AllA$graphit, AllA$mean_val, col= factor(AllA$graphit), pch=19, cex = 3)
  arrows(AllA$graphit, AllA$LCLS, AllA$graphit, AllA$UCLS, length=0.15,lwd=4, angle=90, code=3, col= "black")
  eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
  axis(1, at= c(1:4),labels = c("No CpG \n Syn", " CpG \n Syn", "No CpG \n NonSyn", "CpG \n NonSyn"), mgp=c(3, 1.5, 0))
  axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
  mtext('0', side=2, line=1.5, at=0.0001, las=1.1)

  mtext(nice_name, outer=TRUE, adj=0.55, cex=1.7, line=0.01)
 
  
  plot(jitter(AllT$graphit),AllT$Freq+ 0.0001,log='y',col=factor(AllT$graphit),pch=16,main="T->C",xlab = " ", ylab = "Mutation Frequency", yaxt="n", xaxt = "n", ylim=c(0.0001, .5))
  points(AllT$graphit, AllT$mean_val, col= factor(AllT$graphit), pch=19, cex = 3)
  arrows(AllT$graphit, AllT$LCLS, AllT$graphit, AllT$UCLS, length=0.15, lwd = 4, angle=90, code=3, col= "black")
  eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
  axis(1, at= c(1:4),labels = c("No CpG \n Syn", " CpG \n Syn", "No CpG \n NonSyn", "CpG \n NonSyn"),mgp=c(3, 1.5, 0))
  axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
  mtext('0', side=2, line=1.5, at=0.0001, las=1.1)
  
  
  palette(alpha(c("#99FF99","#FF9900"),0.3)) 
  plot(jitter(AllC$graphit, 0.6),AllC$Freq+ 0.0001,log='y',col=factor(AllC$graphit),pch=16,main="C->T",xlab = "Mutation Type", xlim = c(0.7,4.1), ylab = "Mutation Frequency",yaxt="n", xaxt = "n", ylim=c(0.0001, .5))
  
  points(AllC$graphit, AllC$mean_val, col= factor(AllC$graphit), pch=19, cex = 3)
  arrows(AllC$graphit, AllC$LCLS, AllC$graphit, AllC$UCLS, length=0.15,lwd = 4, angle=90, code=3, col= "black")
  eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
  axis(1, at= c(0.9:3.9),labels = c("No CpG \n Syn", " CpG \n Syn", "No CpG \n NonSyn", "CpG \n NonSyn"),mgp=c(3, 1.5, 0))
  axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
  mtext('0', side=2, line=1.5, at=0.0001, las=1.1) 
  
  
  plot(jitter(AllG$graphit, 0.6),AllG$Freq+ 0.0001,log='y',col=factor(AllG$graphit),pch=16,main="G->A",xlab = "Mutation Type", xlim = c(0.7,4.1), ylab = "Mutation Frequency", yaxt="n", xaxt = "n", ylim=c(0.0001, .5))
  points(AllG$graphit, AllG$mean_val, col= factor(AllG$graphit), pch=19, cex = 3)
  arrows(AllG$graphit, AllG$LCLS, AllG$graphit, AllG$UCLS, length=0.15, lwd=4,angle=90, code=3, col= "black")
  eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
  axis(1, at= c(0.9:3.9),labels = c("No CpG \n Syn", " CpG \n Syn", "No CpG \n NonSyn", "CpG \n NonSyn"),mgp=c(3, 1.5, 0))
  axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
  mtext('0', side=2, line=1.5, at=0.0001, las=1.1)
  dev.off()
} 
