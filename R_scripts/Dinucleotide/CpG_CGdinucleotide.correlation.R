library(colorspace)
library(ggplot2)
library(ggthemes)
library(reshape2)
colors<-qualitative_hcl(6, palette="Dark3")

cpgMF<-read.csv("output/consensus/alldatapoints.csv", stringsAsFactors = F, row.names = 1)
dfA<-cpgMF[,c("Virus","AsynNC_C")]
dfT<-cpgMF[,c("Virus","TsynNC_C")]
colnames(dfA)[2]<-"CpG.Cost.A"
colnames(dfT)[2]<-"CpG.Cost.T"

df<-merge(dfA, dfT, by="Virus")

di<-read.csv("output/Dinuc/CG_Summary.csv", stringsAsFactors = F, row.names = 1)
di<-di[,c("Virus", "cg.rho")]

df2<-merge(df, di, by="Virus")


cor.test(log(df2$CpG.Cost.T), df2$cg.rho,method = "spearman")
#S = 17692, p-value = 0.004444
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#    rho 
#-0.4335953 

cor.test(log(df2$CpG.Cost.A), df2$cg.rho,method = "spearman")
#S = 15904, p-value = 0.06403
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#    rho 
#-0.2887124  


df3<-melt(df2)
rho<-df2[,c("Virus","cg.rho")]
rho<-rbind(rho,rho)
n1<-nrow(rho)
df.plot<-cbind(df3[1:n1,], rho[,2])
colnames(df.plot)[4]<-"CG.rho"
df.plot$variable<-factor(df.plot$variable, levels=c("CpG.Cost.A","CpG.Cost.T" ))

cor.test(df.plot$value, df.plot$CG.rho,method = "spearman" )
#S = 135240, p-value = 0.000547
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#    rho 
#-0.3692594

text1<-paste("A:rho == -0.289")
text2<-paste("T:rho == -0.433")
text<-paste("rho = -0.37, P=0.0005")

ggplot(df.plot, aes(x=CG.rho, y=value, color=variable))+
    scale_y_continuous(trans="log10",breaks = c(1,2,5, 10,20, 50, 100, 200),labels=c(1,2,5, 10,20, 50, 100, 200))+
    geom_point(stat="identity")+
    theme(legend.title = element_blank())+
    scale_color_manual(values=colors[c(1,5)], labels=c("A","T"))+
    theme_bw()+
    ylab("Cost of CpG creating mutations")+
    xlab("Dinucleotide over/udner-representation")+
    theme(panel.grid.minor=element_blank(), panel.grid.major = element_blank())+
    geom_smooth(method='lm',formula=y~x, se=F, size=.3, linetype=2)+
    theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(linetype = 0, size=2)))+
    annotate(geom="text",x=0.75, y=200, label=text, hjust=0, size=2)

ggsave("output/Dinuc/CPG_correlation.plot1.pdf", width = 5,height = 4)  


ggplot(df.plot, aes(x=CG.rho, y=value, color=variable))+
    scale_y_continuous(trans="log10",breaks = c(1,2,5, 10,20, 50, 100, 200),labels=c(1,2,5, 10,20, 50, 100, 200))+
    geom_point(stat="identity")+
    theme(legend.title = element_blank())+
    scale_color_manual(values=colors[c(1,5)], labels=c("A","T"))+
    theme_bw()+
    ylab("Cost of CpG creating mutations")+
    xlab("Dinucleotide over/udner-representation")+
    theme(panel.grid.minor=element_blank(), panel.grid.major = element_blank())+
    geom_smooth(method='lm',formula=y~x, se=F, size=.3, linetype=2)+
    theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(linetype = 0, size=2)))+
    annotate(geom="text",x=0.75, y=200, label=text1, hjust=0, size=2, parse=T)+
    annotate(geom="text",x=0.75, y=150, label=text2, hjust=0, size=2, parse=T)
ggsave("output/Dinuc/CPG_correlation.plot2.pdf", width = 5,height = 4)  
