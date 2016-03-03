rxn2sybil<-function(reaction){
  coeficients <- function(met) {
    regmatches(met, gregexpr('^[0-9]+[[:blank:]]', met))
  }

  if(grepl("[[:blank:]]=>[[:blank:]]",reaction)){
    reactant <- sub("(.*) => (.*)","\\1",reaction)
    products <- sub("(.*) => (.*)","\\2",reaction)
    reactant <- strsplit(reactant," + ",fixed = TRUE)[[1]]
    products <- strsplit(products," + ",fixed = TRUE)[[1]]
    reactant <- gsub("\\+","",reactant)
    products <- gsub("\\+","",products)
    r_coefic <- as.numeric(sapply(reactant, coeficients))
    r_coefic[is.na(r_coefic)]<-1
    p_coefic <- as.numeric(sapply(products, coeficients))
    p_coefic[is.na(p_coefic)]<-1
    reactant <- gsub("^[[:digit:]]+[[:blank:]]","",reactant)
    products <- gsub("^[[:digit:]]+[[:blank:]]","",products)
    reactant <- mapply(function(c,m){paste("(",c,") ", m, sep = "")}, c=r_coefic, m=reactant)
    products <- mapply(function(c,m){paste("(",c,") ", m, sep = "")}, c=p_coefic, m=products)
    reaction <- paste(paste0(reactant,collapse = " + "), paste0(products,collapse = " + "),sep = " --> ")
  }else{
    reactant <- sub("(.*) <=> (.*)","\\1",reaction)
    products <- sub("(.*) <=> (.*)","\\2",reaction)
    reactant <- strsplit(reactant," + ",fixed = TRUE)[[1]]
    products <- strsplit(products," + ",fixed = TRUE)[[1]]
    reactant <- gsub("\\+","",reactant)
    products <- gsub("\\+","",products)
    r_coefic <- as.numeric(sapply(reactant, coeficients))
    r_coefic[is.na(r_coefic)]<-1
    p_coefic <- as.numeric(sapply(products, coeficients))
    p_coefic[is.na(p_coefic)]<-1
    reactant <- gsub("^[[:digit:]]+[[:blank:]]","",reactant)
    products <- gsub("^[[:digit:]]+[[:blank:]]","",products)
    reactant <- mapply(function(c,m){paste("(",c,") ", m, sep = "")}, c=r_coefic, m=reactant)
    products <- mapply(function(c,m){paste("(",c,") ", m, sep = "")}, c=p_coefic, m=products)
    reaction <- paste(paste0(reactant,collapse = " + "), paste0(products,collapse = " + "),sep = " <==> ")
  }
  gsub("\\[s\\]","\\[e\\]",reaction)
}
