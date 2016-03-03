react_compartment <- function(reaction){
  reaction <- rxn2sybil(reaction)
  if(grepl("[[:blank:]]<==>[[:blank:]]",reaction)){
    reaction <- strsplit(reaction," <==> ",fixed = TRUE)[[1]]
  } else {
    reaction <- strsplit(reaction," --> ",fixed = TRUE)[[1]]}
  reaction <- unlist(strsplit(reaction," + ",fixed = TRUE))
  comp<-unlist(regmatches(reaction, gregexpr('[[:punct:]][[:alpha:]]+?\\_?[[:alpha:]]?[[:punct:]]$', reaction)))
  comp<-gsub("\\[","",comp)
  comp<-gsub("\\]","",comp)
  paste(unique(comp),collapse = ", ")
}
