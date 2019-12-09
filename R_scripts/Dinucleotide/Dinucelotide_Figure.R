library(msa)
library(ape)
library(car)
library(seqinr)
library(purrr)
library(dplyr)
library(tidyverse)
library(colorspace)
library(ggplot2)
library(ggthemes)
library(emojifont)
colors<-qualitative_hcl(6, palette="Dark3")

#create summary table of dinucleotide rho values
dint<-list.files("output/ancestor/Dinuc/", pattern = "Dinuc.freq")
Di<-list()
for (i in 1:length(dint)){
    dt<-read.csv(paste0("output/ancestor/Dinuc/", dint[i]), stringsAsFactors = F, row.names = 1)
    dt<-dt[, c("diNT", "Rho")]
    fname<-substr(dint[i], start=12, stop=(nchar(dint[i])-4))
    Di[[i]]<-dt
    names(Di)[i]<-fname
}
for (i in 1:length(Di)) {
    colnames(Di[[i]])<-c("diNt",paste0(names(Di[i])))
}

nt.Rho<-Di%>% purrr::reduce(full_join, by='diNt')    
write.csv(nt.Rho, "output/ancestor/Dinuc/Rho_summary.csv")
#####

#create CpG dinucloetide Rho value summary across all viruses
cpg<-data.frame(virus=colnames(nt.Rho)[2:ncol(nt.Rho)])
cpg$cg.rho<-as.numeric(nt.Rho[which(nt.Rho$diNt=="cg"), 2:ncol(nt.Rho)])
write.csv(cpg, "output/ancestor/Dinuc/CG_rho_values.csv")

# Attach RNA or DNA info. 
info<-read.csv("output/Dinuc/Virus_info.csv", row.names = 1)
cpg<-merge(cpg, info, by="virus")
cpg$Type<-factor(cpg$Type, levels=c("RNA", "DNA"))
cpg<-cpg[order(cpg$Type),]
cpg$Virus<-factor(cpg$Virus, level=paste0(cpg$Virus))
write.csv(cpg, "output/ancestor/Dinuc/CG_summary.csv")

n<-nrow(cpg[cpg$Type=="RNA",])


ggplot(cpg, aes(x=Virus, y=cg.rho, color=Type))+geom_point()+
    geom_point()+
    theme(axis.text.x = element_text(size = 6, angle = 90))+
    scale_color_manual(values=colors[c(1,4)])+
    theme(axis.title.x=element_blank())+ylab("Rho statistic")+
    theme_bw()+
    theme(axis.text.x = element_text(size = 7, angle = 90, hjust=1, color="black"))+
    xlab("")+annotate("rect", xmin=0, xmax=9.5,ymax=Inf,ymin=-Inf, fill="gray90", alpha=0.1 , color=NA)+
    geom_hline(yintercept = 1.23, color = "gray30", size=.5,linetype=1)+
    geom_hline(yintercept = 0.79, color = "gray30", size=.5,linetype=1)+
    geom_vline(xintercept = c(1:(nrow(cpg)-1))+0.5, color = "gray70", size=.2)+
    theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(size = 0.2, colour = "grey30", linetype=2))+  
    annotate("rect", xmin=0, xmax=n+0.5,ymax=Inf,ymin=-Inf, alpha=0.1)+
    geom_text(label=paste0("\u2193 CpG \n  under-represented"), x= nrow(cpg)+1, y=0.75, size=2.5, color="gray20",hjust=0)+
    geom_text(label=paste0("\u2191 CpG \n  over-represented"), x= nrow(cpg)+1, y=1.25, size=2.5, color="gray20",hjust=0)+
    coord_cartesian(clip = 'off') +
    theme(plot.margin = unit(c(.5, 4, .5, .5), "cm"),legend.position = c(1.1, 0.1), legend.title=element_text(size=9),
          legend.text = element_text(size=8), legend.key.size = unit(0.3, "cm"))

#top, right, bottom, and left  
ggsave("output/ancestor/Dinuc/Rho_Figure.pdf", width=9, height=4)







