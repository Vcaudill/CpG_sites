#function for syn/non/nonsence
#This function intakes a dataframe with the columns MUTAA, WTAA_consensus , and TypeOfSite created. Then evaluates the values in the WTAA and MUTAA columns to determine the value for the TypeOfSite column
synFunction <- function(df) {
  for (h in 1:nrow(df)){ #looks at each row in the dataframe df
    if(df[h,"MUTAA"]== df[h,"WTAA_consensus"]){ #if the value in the MUTAA column for the row of interest is equal to the value in the WTAA of the same row
      df[h,"TypeOfSite"] = "syn" #then insert "syn" for the value in the TypeOfSite column for this row
    }
    if(df[h,"MUTAA_ancestor"]== df[h,"WTAA_ancestor"]){ #if the value in the MUTAA column for the row of interest is equal to the value in the WTAA of the same row
      df[h,"ancestor_TypeOfSite"] = "syn" #then insert "syn" for the value in the TypeOfSite column for this row
    }
    if(df[h,"MUTAA"] != df[h,"WTAA_consensus"]){ #if the value in the MUTAA column for the row of interest is NOT equal to the value in the WTAA of the same row
      if(df[h,"MUTAA"]=="*"){ #and it the MUTAA value equal "*"
        df[h,"TypeOfSite"] = "nonsense"} #then insert "nonsense" for the value in the TypeOfSite column for this row
      else {
        df[h,"TypeOfSite"] = "nonsyn" #else insert "nonsyn" for the value in the TypeOfSite column for this row
      }
    }
    if(df[h,"MUTAA_ancestor"] != df[h,"WTAA_ancestor"]){ #if the value in the MUTAA column for the row of interest is NOT equal to the value in the WTAA of the same row
        if(df[h,"MUTAA_ancestor"]=="*"){ #and it the MUTAA value equal "*"
            df[h,"ancestor_TypeOfSite"] = "nonsense"} #then insert "nonsense" for the value in the TypeOfSite column for this row
        else {
            df[h,"ancestor_TypeOfSite"] = "nonsyn" #else insert "nonsyn" for the value in the TypeOfSite column for this row
        }
     }
  }
return(df)
}



#RM