metabolites<-function(reaction){
  if(grepl("[[:blank:]]-->[[:blank:]]",reaction)){
    reaction <- strsplit(reaction," --> ",fixed = TRUE)[[1]]
  } else {
    reaction <- strsplit(reaction," <==> ",fixed = TRUE)[[1]]}
  reaction <- unlist(strsplit(reaction," + ",fixed = TRUE))
  gsub("\\([[:digit:]]+\\)[[:blank:]]","",reaction)
}
