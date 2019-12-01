
Freq <- function(fasta_file, ancestor_seq){
  library(seqinr) 
  virus_basic <- read.fasta(fasta_file)
  number_of_seqs <- length(virus_basic)
  cat("number_of_seqs:", number_of_seqs, "\n")

  virus_align <- read.alignment(fasta_file, format = "fasta", forceToLower = T)
  virus_consensus <- seqinr :: consensus(virus_align, method = "majority")
  virus_ancestor_seq <- tolower(unlist(strsplit(as.character(ancestor_seq), "")))
  #virus_ancesertal=specfic line in a csv file from victoria
  virus_consensus_matrix <- seqinr :: consensus(virus_align, method = "profile")
  virus_ancestor_matrix <- seqinr :: consensus(virus_align, method = "profile")
  consensus_length <- length(virus_consensus)
  number_column <- seq(1, consensus_length)
  virus_DF <- data.frame("num" = number_column, "Freq" = 0, "wtnt_consensus" = virus_consensus,"aFreq" = 0, "ancestor"= "wtnt_ancestor")
 
  for(x in 1:consensus_length){
    current_base <- virus_consensus[x]
    current_ancestor_base <- virus_ancestor_seq[x]
    current_matrix_base_count <- virus_consensus_matrix[,x]
    current_ancestor_matrix_base_count <- virus_ancestor_matrix[,x]
    ts_count <- 0
    as_count <- 0
    if(current_base == "a"){
      ts_count <- current_matrix_base_count[["g"]]
    }
    if(current_base == "g"){
      ts_count <- current_matrix_base_count[["a"]]
    }
    if(current_base == "c"){
      ts_count <- current_matrix_base_count[["t"]]
    }
    if(current_base == "t"){
      ts_count <- current_matrix_base_count[["c"]]
    }
    if(current_ancestor_base == "a"){
        as_count <- current_ancestor_matrix_base_count[["g"]]
    }
    if(current_ancestor_base == "g"){
        as_count <- current_ancestor_matrix_base_count[["a"]]
    }
    if(current_ancestor_base == "c"){
        as_count <- current_ancestor_matrix_base_count[["t"]]
    }
    if(current_ancestor_base == "t"){
        as_count <- current_ancestor_matrix_base_count[["c"]]
    }
    virus_DF[x, 2] <- ts_count/number_of_seqs
    virus_DF[x, 4] <- as_count/number_of_seqs
    virus_DF$ancestor[x]<-as.character(virus_ancestor_seq[x])
  }
  virus_DF$wtnt_consensus<-as.character(virus_DF$wtnt_consensus)
  
  return(virus_DF)
}

