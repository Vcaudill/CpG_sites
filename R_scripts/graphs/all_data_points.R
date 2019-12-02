# graph for data vs significant 

library(scales)
library(plotrix)
library(sfsmisc)
library(gdata)

tally<- read.csv("output/tally.csv")

seq_info<-read.csv("output/Final_CpG_List.csv")

tally <- merge(tally, seq_info ,by="File")

# #Syn CpG vs Non-CpG
#   
tally$SynCpg<-tally$Consensus_a_Syn_CpG_v_NonCpG+tally$Consensus_t_Syn_CpG_v_NonCpG
tally$NonSynCpg<-tally$Consensus_a_NonSyn_CpG_v_NonCpG+tally$Consensus_t_NonSyn_CpG_v_NonCpG
tally$SynNonSyn<-tally$Consensus_a_Syn_v_NonSyn +tally$Consensus_t_Syn_CpG_v_NonCpG
tally$CpG_total<-tally$ca+tally$tg


palette(alpha(c("red","deepskyblue1","green")))

 
png("output/data_points.png", width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2))#, bg = "darkseagreen1"

plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$SynCpg)],pch=c(19,15,17)[as.factor(tally$SynCpg)],cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")


plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$NonSynCpg)],pch=c(19,15,17)[as.factor(tally$NonSynCpg)],cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")

plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$SynNonSyn)],pch=c(19,15,17)[as.factor(tally$SynNonSyn)],cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# of Nucletides")

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(17,15,19))
dev.off()

#### numver of possible cpg sites

palette(alpha(c("red","deepskyblue1","green")))


png("output/data_CpG_sites.png", width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2))#, bg = "darkseagreen1"

plot(tally$Number_of_Sequences,tally$CpG_total,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$SynCpg)],pch=c(19,15,17)[as.factor(tally$SynCpg)],cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of CpG sites")


plot(tally$Number_of_Sequences,tally$CpG_total,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$NonSynCpg)],pch=c(19,15,17)[as.factor(tally$NonSynCpg)],cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of CpG sites")

plot(tally$Number_of_Sequences,tally$CpG_total,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$SynNonSyn)],pch=c(19,15,17)[as.factor(tally$SynNonSyn)],cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# of CpG sites")

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(17,15,19))








#### simple stats

sum(tally$SynCpg)/(2*nrow(tally)) #% both Syn A&T Costly
sum(tally$NonSynCpg)/(2*nrow(tally)) #% both NonSyn A&T Costly
sum(tally$SynNonSyn)/(2*nrow(tally)) #% of Syn more Costly then NonSyn

sum(tally$A..G.Syn..CpG.v.NonCpG.)/nrow(tally) #% Syn A Costly
sum(tally$A..G.NonSyn..CpG.v.NonCpG)/nrow(tally)#% NonSyn A Costly
sum(tally$A..G.Syn.v.NonSyn)/nrow(tally) #% A Syn more Costly then NonSyn

sum(tally$T..C.Syn..CpG.v.NonCpG.)/nrow(tally)#% Syn T Costly
sum(tally$T..C.NonSyn..CpG.v.NonCpG)/nrow(tally) #% NonSyn T Costly
sum(tally$T..C.Syn.v.NonSyn)/nrow(tally) #% T Syn more Costly then NonSyn




##### other types of logged graphs


 png("output/amount/Syn_5_18.png", width = 6.75, height = 6.75, units = "in", res= 300)
  #par(mfrow=c(2,2))#, bg = "darkseagreen1"
  plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,log='xy',col=factor(tally$SynCpg),pch=19,cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")
  legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
         col= c("deepskyblue1","green","red"), horiz= FALSE, 
         #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
         cex=1, pch = c(16,16,16))
  
  png("output/amount/NonSyn_5_18.png", width = 6.75, height = 6.75, units = "in", res= 300)
  plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,log='xy',col=factor(tally$NonSynCpg),pch=19,cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")
  legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
         col= c("deepskyblue1","green","red"), horiz= FALSE, 
         #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
         cex=1, pch = c(16,16,16))
  png("output/amount/SynNonSyn_5_18.png", width = 6.75, height = 6.75, units = "in", res= 300)
  plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,log='xy',col=factor(tally$SynNonSyn),pch=19,cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# of Nucletides")
  legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
         col= c("deepskyblue1","green","red"), horiz= FALSE, 
         #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
         cex=1, pch = c(16,16,16))
dev.off()

png("output/amount/alllogy_5_18.png", width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2))#, bg = "darkseagreen1"

plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,log='y',col=factor(tally$SynCpg),pch=19,cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(16,16,16))


plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,log='y',col=factor(tally$NonSynCpg),pch=19,cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(16,16,16))

plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,log='y',col=factor(tally$SynNonSyn),pch=19,cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# of Nucletides")
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(16,16,16))
dev.off()


png("output/amount/allnolog_5_18.png", width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2))#, bg = "darkseagreen1"

plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,col=factor(tally$SynCpg),pch=19,cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(16,16,16))


plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,col=factor(tally$NonSynCpg),pch=19,cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(16,16,16))

plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,col=factor(tally$SynNonSyn),pch=19,cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# of Nucletides")
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(16,16,16))
dev.off()

