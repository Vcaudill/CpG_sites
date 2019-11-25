###Get Wildtype amino acid###

getWTAA<-function(df){
  check.integer <- function(N){
    !grepl("[^[:digit:]]", format(N,  digits = 20, scientific = FALSE))
  }
  if (check.integer(nrow(df)/3) == FALSE){
    df <- df[-nrow(df),]
    if (check.integer(nrow(df)/3) == FALSE){
      df <- df[-nrow(df),]
    }
    
  }
  #Assign consensus to a variable
  cons =  df$wtnt
  
  #Create empty vector for wildtype amino acid
  WTAA_consensus <- c()
  
  #Loop for translating consensus
  for(x in seq(1, length(cons) - 2, 3)){
    codon <- c(cons[x], cons[x+1], cons[x+2])
    new_AA <- seqinr::translate(codon)
    WTAA_consensus [x] <- new_AA
    WTAA_consensus [x+1] <- new_AA
    WTAA_consensus [x+2] <- new_AA
  }
  
  #Create "WTAA" column if not already
  if (length(which(names(df)=="WTAA_consensus "))==0){
    df$WTAA_consensus =0}
  
  #Insert value into column
  df$WTAA_consensus <-WTAA_consensus 
  
  #Return the data frame
  return(df)
}
