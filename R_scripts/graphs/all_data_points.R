# graph for data vs significant 

library(scales)
library(plotrix)
library(sfsmisc)
library(gdata)
# Virus_info<- read.csv("data/CpG_List.csv")
# 
# data_points_df<- read.csv("output/alldatapoints.csv")
# #data_points_df=data_points_df[-(27), ] 
# data_points_df <- drop.levels(data_points_df[-(27), ]) 
# tally<- read.csv("data/recount2.csv")
tally<- read.csv("output/All_Data/Data_wilcox_tally_results/CpG_Tally_923.csv")

# data_points_df$SynCpg=0
# data_points_df$NonSynCpg=0
# data_points_df$SynNonSyn=0
# for (i in 1:nrow(data_points_df)) {
#  treat=data_points_df$Virus[i]
#   
#   print(treat)
#  # treat=factor(data_points_df$Virus[i])
#   for (j in 1:nrow(tally)){
#     # splitnameAll<-unlist(strsplit(as.character(Virus_info$name[j]),".fasta"))
#     # splitname<-unlist(strsplit(as.character(splitnameAll[1]),"_"))
#    truename<-tally$Virus[j]
#     print(truename)
#    # truename=factor(tally$Virus[i])
#     if(treat==truename){
#       data_points_df$SynCpg[i]<-tally$A..G.Syn..CpG.v.NonCpG.[j]+tally$T..C.Syn..CpG.v.NonCpG.[j]
#       data_points_df$NonSynCpg[i]<-tally$A..G.NonSyn..CpG.v.NonCpG[j]+tally$T..C.NonSyn..CpG.v.NonCpG[j]
#       data_points_df$SynNonSyn[i]<-tally$A..G.Syn.v.NonSyn[j]+tally$T..C.Syn.v.NonSyn[j]
#       
#       
#     }
#   } 
# }
# 
# #Syn CpG vs Non-CpG
#   
tally$SynCpg<-tally$A..G.Syn..CpG.v.NonCpG.+tally$T..C.Syn..CpG.v.NonCpG.
tally$NonSynCpg<-tally$A..G.NonSyn..CpG.v.NonCpG+tally$T..C.NonSyn..CpG.v.NonCpG
tally$SynNonSyn<-tally$A..G.Syn.v.NonSyn+tally$T..C.Syn.v.NonSyn

palette(alpha(c("red","deepskyblue1","green")))

 
png("output/All_Data/Num_of_data_points/alllogxy_9_25.png", width = 6.75, height = 6.75, units = "in", res= 300)
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

