met_compartment<-function(m_name,metabolites){
  m_name <- gsub("\\]","*",m_name)
  m_name <- gsub("\\[","*",m_name)
  m_name <- gsub("\\(","*",m_name)
  m_name <- gsub("\\)","*",m_name)
  m_name <- gsub("\\*","[[:punct:]]",m_name)
  m_name <- paste0("^",m_name,"[[:punct:]]")
  mets <- metabolites[grep(m_name,metabolites)]
  comp<-unlist(regmatches(mets, gregexpr('\\[[[:alpha:]]\\]$', mets)))
  comp<-gsub("\\[","",comp)
  comp<-gsub("\\]","",comp)
  paste0(unique(comp),collapse = ", ")
}
