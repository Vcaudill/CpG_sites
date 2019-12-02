######## Investigating makesCpG sites: 

## Note : this function assumes that wtnt column is present in the 
##        dataframe and your dataframe that is set as df

# Runs makesCpG_site function which takes in the arguement/input dataframe
CPG_site<-function(df){
  
  # Collapses wtnt_consensus vector sequence into a character string and sets to variable STRING
  paste(df$wtnt_consensus, collapse = '') -> STRING
  paste(df$ancestor, collapse = '') -> ASTRING
  STRING
  
  #Looks for pattern tg in the data STRING, return location of TG sites within the STRING and store in variable TG
  gregexpr(pattern ='tg',STRING ) -> TG
  gregexpr(pattern ='tg',ASTRING ) -> ATG
  TG
  
  # Since TG is a list, insert 'list' TG into a data frame BELL
  BELL <- data.frame(matrix(unlist(TG)))
  ABELL <- data.frame(matrix(unlist(ATG)))
  BELL
  
  # create new column name makesCpG with values zero
  df$makesCpG<- 0
  df$ancestor_makesCpG<- 0
  
  # Inserting value 1 in the column "makesCpG" and using the values found dataframe BELL as the row #
  df[BELL[,1],"makesCpG"] <- 1
  df[ABELL[,1],"ancestor_makesCpG"] <- 1
  
  
  # For CA sites: 
  
  #Looks for pattern ca in the data STRING, return location of CA sites within the STRING and store in variable CA
  gregexpr(pattern ='ca',STRING ) -> CA
  gregexpr(pattern ='ca',ASTRING ) -> ACA
  CA
  
  # Since TG is a list, insert 'list' TG into a data frame BELL
  CASITES <- data.frame(matrix(unlist(CA)))
  ACASITES <- data.frame(matrix(unlist(ACA)))
  CASITES
  
  # Inserting value 1 in the column "makesCpG" and using the values found dataframe BELL as the row #
  df[CASITES[,1]+1,"makesCpG"] <- 1
  df[ACASITES[,1]+1,"ancestor_makesCpG"] <- 1
  #returns new dataframe
  return(df)
}

# Calls the makesCpG_site function with the dataframe as the argument Old dataframe becomes the new dataframe


