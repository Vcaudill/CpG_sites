
library(scales)
library(plotrix)
library(sfsmisc)
library(varhandle)
#makesCpG<->ancestor_makesCpG
#wtnt_consensus<->ancestor
#TypeOfSite<->ancestor_TypeOfSite
mean_or_median<- "mean"
#Function to help create errorbars
sem<-function(x){
    return(sd(x,na.rm = FALSE)/sqrt(length(x)))
}
#cycle through csv

my.list<- list.files("output/Csv", pattern=".csv")
data_points = data.frame("order"= 1:length(my.list), "Virus"= 0,'total'=0,"AsynNC_C" =0,'AsynC_LCLS' = 0, 'AsynNC_LCLS'= 0, 'AsynC_UCLS'= 0, 'AsynNC_UCLS'= 0,'TsynNC_C'=0,'TsynC_LCLS' = 0, 'TsynNC_LCLS' = 0,  'TsynC_UCLS'= 0,  'TsynNC_UCLS' = 0 )
Virus_info<-read.csv("output/Final_CpG_List.csv")
count = 1
for (i in 1:length(my.list)){  
    data<-read.csv(paste0("output/Csv/",my.list[i]))
    Name<-Virus_info[grep(strsplit(my.list[i],".csv"), Virus_info$File), ]
    cpg.y<-subset(data, ancestor_makesCpG==1)
    cpg.n<-subset(data, ancestor_makesCpG==0)
    #subset further into letters nuclotideCpgforming or nucotideNonGpg
    AC<-subset(cpg.y, ancestor=='a')
    ANC<-subset(cpg.n, ancestor=='a')
    TC<-subset(cpg.y, ancestor=='t')
    TNC<-subset(cpg.n, ancestor=='t')
    
    
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
    for (j in 1:length(AllA$ancestor_makesCpG)) {
        if (AllA$ancestor_makesCpG[j] == 1 && AllA$ancestor_TypeOfSite[j] == "syn") {
            AllA_mean_value_syn_CpG <- stats(AllA$aFreq[(which(AllA$ancestor_makesCpG == 1 & AllA$ancestor_TypeOfSite == "syn") )])
            AllA_sem_vals_syn_CpG<-error_bar(AllA$aFreq[(which(AllA$ancestor_makesCpG == 1 & AllA$ancestor_TypeOfSite == "syn") )])
            
        }
        
        if (AllA$ancestor_makesCpG[j] == 0 && AllA$ancestor_TypeOfSite[j] == "syn") {
            AllA_mean_value_syn_nCpG <- stats(AllA$aFreq[(which(AllA$ancestor_makesCpG == 0 & AllA$ancestor_TypeOfSite == "syn") )])
            
            AllA_sem_vals_syn_nCpG<-error_bar(AllA$aFreq[(which(AllA$ancestor_makesCpG == 0 & AllA$ancestor_TypeOfSite == "syn") )])
        }
       
    }
    
    for (k in 1:length(AllT$ancestor_makesCpG)) {
        if (AllT$ancestor_makesCpG[k] == 1 && AllT$ancestor_TypeOfSite[k] == "syn") {
            AllT_mean_value_syn_CpG <- stats(AllT$aFreq[(which(AllT$ancestor_makesCpG == 1 & AllT$ancestor_TypeOfSite == "syn") )])
            AllT_sem_vals_syn_CpG<-error_bar(AllT$aFreq[(which(AllT$ancestor_makesCpG == 1 & AllT$ancestor_TypeOfSite == "syn") )])
        }
        
        if (AllT$ancestor_makesCpG[k] == 0 && AllT$ancestor_TypeOfSite[k] == "syn") {
            AllT_mean_value_syn_nCpG<- stats(AllT$aFreq[(which(AllT$ancestor_makesCpG == 0 & AllT$ancestor_TypeOfSite == "syn") )])
            AllT_sem_vals_syn_nCpG<-error_bar(AllT$aFreq[(which(AllT$ancestor_makesCpG == 0 & AllT$ancestor_TypeOfSite == "syn") )])
        }
        
    }
    
    # There are the upper and lower limits of the error bar
    # AllA_LCLS = AllA$mean_value - AllA$sem_vals
    # AllA_UCLS = AllA$mean_value + AllA$sem_vals
    # 
    # AllT_LCLS = AllT$mean_value - AllT$sem_vals
    # AllT_UCLS = AllT$mean_value + AllT$sem_vals
    # 
    
    data_points$Virus[count]<- unfactor(Name$name)
    data_points$total[count]<- Name$Number_of_Sequences*Name$Number_of_Nucleotides
    data_points$order[count]<-Name$graph_order
    data_points$AsynNC_C[count]= AllA_mean_value_syn_nCpG/AllA_mean_value_syn_CpG

    data_points$TsynNC_C[count] = AllT_mean_value_syn_nCpG/AllT_mean_value_syn_CpG

    
    data_points$AsynC_LCLS[count]= AllA_mean_value_syn_CpG - AllA_sem_vals_syn_CpG
    data_points$AsynNC_LCLS[count]= AllA_mean_value_syn_nCpG - AllA_sem_vals_syn_nCpG
    data_points$TsynC_LCLS[count] = AllT_mean_value_syn_CpG - AllT_sem_vals_syn_CpG
    data_points$TsynNC_LCLS[count] = AllT_mean_value_syn_nCpG - AllT_sem_vals_syn_nCpG 

    data_points$AsynC_UCLS[count]= AllA_mean_value_syn_CpG + AllA_sem_vals_syn_CpG
   
    data_points$AsynNC_UCLS[count]= AllA_mean_value_syn_nCpG + AllA_sem_vals_syn_nCpG

    data_points$TsynC_UCLS[count] = AllT_mean_value_syn_CpG + AllT_sem_vals_syn_CpG
    
    data_points$TsynNC_UCLS[count] = AllT_mean_value_syn_nCpG + AllT_sem_vals_syn_nCpG 
    
    count = count +1
}


df <- apply(data_points,2,as.character)
write.csv(df, file = "output/supplenentary_based_on_ancestor/alldatapoints_ancestor.csv")


# graphing 
png("output/supplenentary_based_on_ancestor/Costly_Graph_ancestor_AllR_12_11_2019.png", width = 15, height = 8, units = "in", res= 500)
#--------------------
par(mar=c(0,2,3,2), oma=c(6,4,1,1), mfrow=c(2,1))#, bg = "darkseagreen1"
#changed mar(0,2,3,2) oma(6,4,1,1)
x <- data_points$order 
y <-  data_points$AsynNC_C
plot(x,y, type = "n", log ='y' ,main="Cost of CpG-Creating Mutations (For Ancestor)", xlab=" ", yaxt = "n", ylab="Costly", xaxt = "n", ylim=c(0.02, 300), xlim=c(2, length(my.list) +5.5), las= 1, cex.main=3) 

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
rect(31.5, .0001, 38.5, 1570, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(39.5, .0001, 40.5, 1570, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(42.5, .0001, 50.5, 1570, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))

# aty <- axTicks(2)
# aty <- axTicks(2)
# labels <- sapply(aty,function(i)
#   as.expression(bquote(10^ .(i)))
# )
# axis(2,at=aty,labels=labels)
points(data_points$order-.3, data_points$AsynNC_C, col='red', pch= 19)
#points(data_points$Count -.1, data_points$AnonsynNC_C, col= "green", pch=19)
points(data_points$order+.1, data_points$TsynNC_C, col= "blue", pch=19)
#points(data_points$Count + .3, data_points$TnonsynNC_C, col= "purple", pch=19)
# par(new=TRUE)

# hack: we draw arrows but with very special "arrowheads"
arrows(data_points$order-.3, data_points$AsynNC_LCLS/data_points$AsynC_LCLS, data_points$order -.3, data_points$AsynNC_UCLS/data_points$AsynC_UCLS, length=0.05, angle=90, code=3, col= "red")
#arrows(data_points$Count-.1, data_points$AnonsynNC_LCLS/data_points$AnonsynC_LCLS, data_points$Count -.1, data_points$AnonsynNC_UCLS/data_points$AnonsynC_UCLS, length=0.05, angle=90, code=3, col= "green")
arrows(data_points$order+.1, data_points$TsynNC_LCLS/data_points$TsynC_LCLS, data_points$order+.1, data_points$TsynNC_UCLS/data_points$TsynC_UCLS, length=0.05, angle=90, code=3, col= "blue")
#arrows(data_points$Count+.3, data_points$TnonsynNC_LCLS/data_points$TnonsynC_LCLS, data_points$Count +.3, data_points$TnonsynNC_UCLS/data_points$TnonsynC_UCLS, length=0.05, angle=90, code=3, col= "purple")

axis(2, at = c(0.02, 0.5,1,2,5,10,20,50,100,200), labels = c(0.02, 0.5,1,2,5,10,20,50,100,200),  las=2)
#axis.break(2, 0.007,breakcol="black",style="slash")
#mtext('No nonCpG \n  or  CpG     \n    mutations  ', side=2, line=.005, at=0.002, las=1.1, cex = .7)
#mtext('No nonCpG \n   mutations   ', side=2, line=.005, at=0.0005, las=1.1, cex = .7)

#axis.break(2, 200,breakcol="black",style="slash")
#mtext('2.47e+13', side=2, line=.005, at=300, las=1.1, cex = .9)

mtext('Relative Cost of CpG \n Creating Mutations', side=2, line=3, at=5, las=0, cex = 1.3)

#mtext('No CpG \n mutations ', side=2, line=.005, at=700, las=1.1, cex = .7)

abline(h=c(0.002, 0.0005,0.2, .01,0.5,1,2,5,10,20,50,100, 200), col="grey", lty=c(2,2))

abline(v=c(1.5,2.5,3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,12.5,13.5,14.5,16.5,17.5,18.5,19.5,15.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5,41.5,42.5), col="grey", lty=c(1))
#abline(v=c(6.5,10.5,14.5,20.5, 27.5, 41.5), col = "darkgreen", lwd=2)
abline(v=c(4.5,6.5,7.5,11.5,17.5,24.5,25.5,26.5,28.5,29.5,30.5,31.5,38.5,39.5,40.5,42.5), col = "bisque4", lwd=2)


abline(h = 1, col ="darkslategrey", lwd = 2)
#rect(0,800, 800,800, col = "darkseagreen1")
# xlab="Virus "
text(length(my.list)+5, 3, " CpG Mutation \n More Costly", cex = 1, font = 2)
text(length(my.list)+2, 3, "↑", cex = 3, font = 2)
text(length(my.list)+5, .4, " CpG Mutation \n Less Costly", cex = 1, font = 2)
text(length(my.list)+2, .4, "↓", cex = 3, font = 2)
#axis(1, at=1:length(my.list), labels=data_points$Virus, las= 2)
legend('topright', legend=c("A -> G Syn","", "T -> C Syn", " "),
       col=c("red", "white", "blue", "white"), lty=1, lwd= 3, cex = 1, pt.cex = 999)


par(mar=c(10,2,0,2))
plot(data_points$order,data_points$total, ylim=c(60000, 45000000), type = "n", log ='y' , xlab=" ", yaxt = "n",  xaxt = "n", xlim=c(2, length(my.list) +5.5), las= 1, cex.main=3) 
mtext('Amount of \n Data', side=2, line=3, at=1100000, las=0, cex = 1.3)

abline(v=c(1.5,2.5,3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,12.5,13.5,14.5,16.5,17.5,18.5,19.5,15.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5), col="grey", lty=c(1))
abline(v=c(4.5,6.5,7.5,11.5,17.5,24.5,25.5,26.5,28.5,29.5,30.5,31.5,38.5,39.5,40.5,42.5), col = "bisque4", lwd=2)

abline(h = 55988780, col ="darkslategrey", lwd = 2, lty = 2)
axis(2, at = c(150000,1000000,5000000,25000000), labels = c("150k","1 mil","5 mil","25 mil"),  las=2)
abline(h=c(150000,1000000,5000000,25000000), col="grey", lty=c(2,2))

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
rect(31.5, .0001, 38.5, 558988780, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(39.5, .0001, 40.5, 558988780, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(42.5, .0001, 50.5, 558988780, density = NULL, angle = 45,
     col = rgb(211/255,211/255,211/255, alpha=.3), border = NULL, lty = par("lty"), lwd = par("lwd"))
points(data_points$order,data_points$total, pch = 16)

axis(1, at=data_points$order, labels=data_points$Virus, las= 2)

dev.off()
# add horsontal gray lines

