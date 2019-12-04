###Get mutated amino acid###
getMUTAA <- function(df){
  
  #Assign consensus to a variable
  cons = df$wtnt_consensus
  acons = df$ancestor
  
  #Create empty vector for mutated amino acid
  MUTAA <- c()
  MUTAA_ancestor <- c()
  
  #Loop for mutated codon
  for(x in seq(1, length(cons), 3)){
    codon <- c(cons[x], cons[x+1], cons[x+2])
    mutated_codon <- codon
    if(codon[1] == "a"){
      mutated_codon <- replace(x=mutated_codon, values=c("g", codon[2], codon[3]))
    }
    if(codon[1] == "g"){
      mutated_codon <- replace(x=mutated_codon, values=c("a", codon[2], codon[3]))
    }
    if(codon[1] == "c"){
      mutated_codon <- replace(x=mutated_codon, values=c("t", codon[2], codon[3]))
    }
    if(codon[1] == "t"){
      mutated_codon <- replace(x=mutated_codon, values=c("c", codon[2], codon[3]))
    }
    MUTAA[x] <- translate(mutated_codon)

    if(codon[2] == "a"){
      mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "g", codon[3]))
    }
    if(codon[2] == "g"){
      mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "a", codon[3]))
    }
    if(codon[2] == "c"){
      mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "t", codon[3]))
    }
    if(codon[2] == "t"){
      mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "c", codon[3]))
    }
    MUTAA[x+1] <- translate(mutated_codon)
    
    if(codon[3] == "a"){
      mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "g"))
    }
    if(codon[3] == "g"){
      mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "a"))
    }
    if(codon[3] == "c"){
      mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "t"))
    }
    if(codon[3] == "t"){
      mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "c"))
    }
    MUTAA[x+2] <- translate(mutated_codon)
  }
  
  for(i in seq(1, length(acons), 3)){   
    #ancestor
    acodon <- c(acons[i], acons[i+1], acons[i+2])
    mutated_a_codon <- acodon
    if(acodon[1] == "a"){
        mutated_a_codon <- replace(x=mutated_a_codon, values=c("g", acodon[2], acodon[3]))
    }
    if(acodon[1] == "g"){
        mutated_a_codon <- replace(x=mutated_a_codon, values=c("a", acodon[2], acodon[3]))
    }
    if(acodon[1] == "c"){
        mutated_a_codon <- replace(x=mutated_a_codon, values=c("t", acodon[2], acodon[3]))
    }
    if(acodon[1] == "t"){
        mutated_a_codon <- replace(x=mutated_a_codon, values=c("c", acodon[2], acodon[3]))
    }
    MUTAA_ancestor[i] <- translate(mutated_a_codon)
    
    if(acodon[2] == "a"){
        mutated_a_codon <- replace(x=mutated_a_codon, values=c(acodon[1], "g", acodon[3]))
    }
    if(acodon[2] == "g"){
        mutated_a_codon <- replace(x=mutated_a_codon, values=c(acodon[1], "a", acodon[3]))
    }
    if(acodon[2] == "c"){
        mutated_a_codon <- replace(x=mutated_a_codon, values=c(acodon[1], "t", acodon[3]))
    }
    if(acodon[2] == "t"){
        mutated_a_codon <- replace(x=mutated_a_codon, values=c(acodon[1], "c", acodon[3]))
    }
    MUTAA_ancestor[i+1] <- translate(mutated_a_codon)
    
    if(acodon[3] == "a"){
        mutated_a_codon <- replace(x=mutated_a_codon, values=c(acodon[1], acodon[2], "g"))
    }
    if(acodon[3] == "g"){
        mutated_a_codon <- replace(x=mutated_a_codon, values=c(acodon[1], acodon[2], "a"))
    }
    if(acodon[3] == "c"){
        mutated_a_codon <- replace(x=mutated_a_codon, values=c(acodon[1], acodon[2], "t"))
    }
    if(acodon[3] == "t"){
        mutated_a_codon <- replace(x=mutated_a_codon, values=c(acodon[1], acodon[2], "c"))
    }
    MUTAA_ancestor[i+2] <- translate(mutated_a_codon)
  }
  
  #Create "MUTAA" category if not already
  if (length(which(names(df)=="MUTAA"))==0){
    df$MUTAA=0}
  if (length(which(names(df)=="MUTAA_ancestor"))==0){
      df$MUTAA_ancestor=0}
  
  #Insert value into column
  df$MUTAA<-MUTAA
  df$MUTAA_ancestor<-MUTAA_ancestor
  
  #Return data frame
  return(df)
}