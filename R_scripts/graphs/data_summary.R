# graph for data vs significant 

library(scales)
library(plotrix)
library(sfsmisc)
library(gdata)

tally<- read.csv("output/Tally.csv")

extra_ifo<-read.csv("output/Extra_CpG_info.csv")
extra_ifo$Count_of_CpG_mutations<-0
extra_ifo$ancestor_Count_of_CpG_mutations<-0

my.list<- list.files("output/Csv", pattern=".csv")

for (i in 1:length(my.list)){  
    data<-read.csv(paste0("output/Csv/",my.list[i]))
    Name<-extra_ifo[grep(strsplit(my.list[i],".csv"), extra_ifo$File), ]
    CpG_mutation_data<-subset(subset(subset(data, makesCpG==1), TypeOfSite == "syn"), Freq>0)
    if (Name$File == extra_ifo$File[i]){
        extra_ifo$Count_of_CpG_mutations[i]<-nrow(CpG_mutation_data)}
    
    data_a<-read.csv(paste0("output/Csv/",my.list[i]))
    Name_a<-extra_ifo[grep(strsplit(my.list[i],".csv"), extra_ifo$File), ]
    a_CpG_mutation_data<-subset(subset(subset(data_a, ancestor_makesCpG==1), ancestor_TypeOfSite == "syn"), aFreq>0)
    if (Name_a$File == extra_ifo$File[i]){
    extra_ifo$ancestor_Count_of_CpG_mutations[i]<-nrow(a_CpG_mutation_data)}
    }

tally <- merge(extra_ifo, tally ,by="File")
# #Syn CpG vs Non-CpG
# consensus
tally$SynCpg<-tally$Consensus_a_Syn_CpG_v_NonCpG+tally$Consensus_t_Syn_CpG_v_NonCpG
tally$NonSynCpg<-tally$Consensus_a_NonSyn_CpG_v_NonCpG+tally$Consensus_t_NonSyn_CpG_v_NonCpG
tally$SynNonSyn<-tally$Consensus_a_Syn_v_NonSyn +tally$Consensus_t_Syn_CpG_v_NonCpG
tally$CpG_total<-tally$ca+tally$tg
# ancestor
tally$ASynCpg<-tally$Ancestor_a_Syn_CpG_v_NonCpG+tally$Ancestor_t_Syn_CpG_v_NonCpG
tally$ANonSynCpg<-tally$Ancestor_a_NonSyn_CpG_v_NonCpG+tally$Ancestor_t_NonSyn_CpG_v_NonCpG
tally$ASynNonSyn<-tally$Ancestor_a_Syn_v_NonSyn +tally$Ancestor_t_Syn_CpG_v_NonCpG
tally$A_CpG_total<-tally$ancestor_ca +tally$ancestor_tg

palette(alpha(c("red","deepskyblue1","green")))

png("output/data_summary/data_points.png", width = 6.75, height = 6.75, units = "in", res= 300)
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

#### number of possible cpg sites

palette(alpha(c("red","deepskyblue1","green")))


png("output/data_summary/data_CpG_sites.png", width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2))#, bg = "darkseagreen1"

plot(tally$Number_of_Sequences,tally$CpG_total,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$SynCpg)],pch=c(19,15,17)[as.factor(tally$SynCpg)],cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of possible CpG sites")


plot(tally$Number_of_Sequences,tally$CpG_total,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$NonSynCpg)],pch=c(19,15,17)[as.factor(tally$NonSynCpg)],cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of possible CpG sites")

plot(tally$Number_of_Sequences,tally$CpG_total,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$SynNonSyn)],pch=c(19,15,17)[as.factor(tally$SynNonSyn)],cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# of possible CpG sites")

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(17,15,19))

dev.off()


#### number of cpg mutations

palette(alpha(c("red","deepskyblue1","green")))


png("output/data_summary/data_CpG_mutations.png", width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2))#, bg = "darkseagreen1"

plot(tally$Number_of_Sequences,tally$Count_of_CpG_mutations,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$SynCpg)],pch=c(19,15,17)[as.factor(tally$SynCpg)],cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# CpG mutations")


plot(tally$Number_of_Sequences,tally$Count_of_CpG_mutations,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$NonSynCpg)],pch=c(19,15,17)[as.factor(tally$NonSynCpg)],cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# CpG mutations")

plot(tally$Number_of_Sequences,tally$Count_of_CpG_mutations,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$SynNonSyn)],pch=c(19,15,17)[as.factor(tally$SynNonSyn)],cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# CpG mutations")

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(17,15,19))

dev.off()



######################################################################
############################### ancestor ###################################
palette(alpha(c("red","deepskyblue1","green")))

png("output/supplenentary_based_on_ancestor/data_summary/data_points.png", width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2))#, bg = "darkseagreen1"

plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$ASynCpg)],pch=c(19,15,17)[as.factor(tally$ASynCpg)],cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")


plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$ANonSynCpg)],pch=c(19,15,17)[as.factor(tally$ANonSynCpg)],cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of Nucletides")

plot(tally$Number_of_Sequences,tally$Number_of_Nucleotides,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$ASynNonSyn)],pch=c(19,15,17)[as.factor(tally$ASynNonSyn)],cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# of Nucletides")

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(17,15,19))
dev.off()

#### numver of possible cpg sites

palette(alpha(c("red","deepskyblue1","green")))


png("output/supplenentary_based_on_ancestor/data_summary/data_CpG_sites.png", width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2))#, bg = "darkseagreen1"

plot(tally$Number_of_Sequences,tally$A_CpG_total,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$ASynCpg)],pch=c(19,15,17)[as.factor(tally$ASynCpg)],cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of possible CpG sites")


plot(tally$Number_of_Sequences,tally$A_CpG_total,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$ANonSynCpg)],pch=c(19,15,17)[as.factor(tally$ANonSynCpg)],cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of possible CpG sites")

plot(tally$Number_of_Sequences,tally$A_CpG_total,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$ASynNonSyn)],pch=c(19,15,17)[as.factor(tally$ASynNonSyn)],cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# of possible CpG sites")

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("bottomright",legend=c("Significant","Partially Significant","No Significance"),
       col= c("deepskyblue1","green","red"), horiz= FALSE, 
       #lty is type of line used,cex is size of legend, xpd allows legend to lie outside the plot,bty is type of box around legend       
       cex=1, pch = c(17,15,19))

dev.off()

#### numver of  cpg mutations

palette(alpha(c("red","deepskyblue1","green")))


png("output/supplenentary_based_on_ancestor/data_summary/data_CpG_mutations.png", width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2))#, bg = "darkseagreen1"

plot(tally$Number_of_Sequences,tally$ancestor_Count_of_CpG_mutations,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$ASynCpg)],pch=c(19,15,17)[as.factor(tally$ASynCpg)],cex=1.5, main="Syn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of CpG Mutations")


plot(tally$Number_of_Sequences,tally$ancestor_Count_of_CpG_mutations,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$ANonSynCpg)],pch=c(19,15,17)[as.factor(tally$ANonSynCpg)],cex=1.5, main="NonSyn CpG vs NonCpG",xlab = "# of Sequences", ylab = "# of CpG Mutations")

plot(tally$Number_of_Sequences,tally$ancestor_Count_of_CpG_mutations,log='xy',col=c("red","green","deepskyblue1")[as.factor(tally$ASynNonSyn)],pch=c(19,15,17)[as.factor(tally$ASynNonSyn)],cex=1.5, main="Syn vs NonSyn",xlab = "# of Sequences", ylab = "# of CpG Mutations")


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

sum(tally$Consensus_a_Syn_CpG_v_NonCpG)/nrow(tally) #% Syn A Costly
sum(tally$Consensus_a_NonSyn_CpG_v_NonCpG)/nrow(tally)#% NonSyn A Costly
sum(tally$Consensus_a_Syn_v_NonSyn)/nrow(tally) #% A Syn more Costly then NonSyn

sum(tally$Consensus_t_Syn_CpG_v_NonCpG)/nrow(tally)#% Syn T Costly
sum(tally$Consensus_t_NonSyn_CpG_v_NonCpG)/nrow(tally) #% NonSyn T Costly
sum(tally$Consensus_t_Syn_v_NonSyn)/nrow(tally) #% T Syn more Costly then NonSyn

#### simple stats ancestor

sum(tally$ASynCpg)/(2*nrow(tally)) #% both Syn A&T Costly
sum(tally$ANonSynCpg)/(2*nrow(tally)) #% both NonSyn A&T Costly
sum(tally$ASynNonSyn)/(2*nrow(tally)) #% of Syn more Costly then NonSyn

sum(tally$Ancestor_a_Syn_CpG_v_NonCpG)/nrow(tally) #% Syn A Costly
sum(tally$Ancestor_a_NonSyn_CpG_v_NonCpG)/nrow(tally)#% NonSyn A Costly
sum(tally$Ancestor_a_Syn_v_NonSyn)/nrow(tally) #% A Syn more Costly then NonSyn

sum(tally$Ancestor_t_Syn_CpG_v_NonCpG)/nrow(tally)#% Syn T Costly
sum(tally$Ancestor_t_NonSyn_CpG_v_NonCpG)/nrow(tally) #% NonSyn T Costly
sum(tally$Ancestor_t_Syn_v_NonSyn)/nrow(tally) #% T Syn more Costly then NonSyn

