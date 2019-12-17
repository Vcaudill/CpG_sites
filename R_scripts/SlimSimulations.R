
if (TRUE) {
  #setwd("~/Desktop/CpG_Project")
  Sums<-rep(0,10000)
  filepath="../data/dataSLIMsimulations/"
  ListOfFiles=paste("../data/dataSLIMsimulations/",list.files(path = "../data/dataSLIMsimulations/"),sep = "")
  
  for (f in ListOfFiles){
    R<-read.table(f)
    for (i in 1:length(R)){
      if (R[i]>0){Sums[i]=Sums[i]+1}}}
  
  mean(as.numeric(Sums[1:5000]))/length(list.files(path = filepath))
  mean(as.numeric(Sums[5001:10000]))/length(list.files(path = filepath))
  mean(as.numeric(Sums[5001:10000]))/mean(as.numeric(Sums[1:5000]))
  
  
  library(scales)
  library(plotrix)
  library(sfsmisc)
  
  sem<-function(x){
    return(sd(x,na.rm = FALSE)/sqrt(length(x)))
  }
}

png("../output/SlimSimulations/SLIMfreq.png", width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(1,1), mar=c(4.1, 4.1, 1.9, 0.8),oma=c(0.1,0.1,1.5,0.1)) 
palette(alpha(c("#99FF99","#9999FF","#FF9900","#FF3300"),0.3))

sem_1 = sd(Sums[1:5000],na.rm = FALSE)/sqrt(length(Sums[1:5000]))
sem_1 = sem(Sums[1:5000])
mean_1 = mean(as.numeric(Sums[1:5000]))/length(list.files(path = filepath))
sem_2 = sd(Sums[5001:10000],na.rm = FALSE)/sqrt(length(Sums[1:5000]))
sem_2 = sem(Sums[5001:10000])
mean_2 = mean(as.numeric(Sums[5001:10000]))/length(list.files(path = filepath))

pos=2
plot(jitter(rep(pos,5000), amount = .3),
     Sums[1:5000]/(length(list.files(path = filepath))) + 0.0001,
                   log='y',col=alpha("#9999FF",0.3),pch=16, main="Simulation results",xlab = " ", ylab = "Mutation Frequency", yaxt="n", xaxt="n", xlim=c(0.5,2.5),ylim=c(0.0001, 1.2))

points(pos, mean_1, col= "#9999FF", pch=19, cex = 3)
arrows(pos, mean_1+sem_1, pos, max(0.0001,mean_1-sem_1), length=0.25, lwd=4, angle=90, code=3, col= "black")
eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
axis(1, at= c(1:2),labels = c("No CpG \n Syn", " CpG \n Syn"), mgp=c(3, 1.5, 0))
axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
mtext('0', side=2, line=1.5, at=0.0001, las=1.1)


pos=1
points(jitter(rep(pos,5000), amount = 0.3),
     Sums[5001:10000]/(length(list.files(path = filepath))) + 0.0001,col=1,pch=16)
points(pos, mean_2, col= "#99FF99", pch=19, cex = 3)
arrows(pos, mean_2+sem_2, pos, max(0.0001,mean_2-sem_2), length=0.25, lwd=4, angle=90, code=3, col= "black")

dev.off()
