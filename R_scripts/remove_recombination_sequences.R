
library(seqinr)
fasta<-"RotavirusA_VP6_alignment.fasta"
fasta_file <-paste("data/",fasta,sep = "")
virus <- read.fasta(fasta_file, as.string = TRUE)
IDs <-names(virus)

ID_list<-"KX268780.1	GQ477091.2	KJ482497.1	HM348746.1	JN014005.1	JN014003.1	HM988972.1	KU714448.1	HQ611009.2	FJ685614.1	KJ753780.1	JX040423.1	KC140589.1	GU199496.1	LC433811.1	LC433800.1	KF740531.1	MG676136.1	KP882441.1"

#adds | in between space
remove_space<-gsub("[[:blank:]]+", "|", ID_list)

for (i in 1:length(virus)){ 
    if(grepl(remove_space,format(IDs[i]))==TRUE){
        next()
    }
    write.fasta(virus[i], open = "a", names = names(virus[i]), file = paste("data/recombination_removed/",fasta,sep=""))
    
}
