reversibility <- function(reaction){
  reaction <- rxn2sybil(reaction)
  if(grepl("[[:blank:]]-->[[:blank:]]",reaction)){
    return("irreversible")
  } else {
    return("reversible")
  }
}
