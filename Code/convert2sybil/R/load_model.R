load_model <- function(file) {
  model<-gdata::read.xls(file,sheet = 1)
  as.data.frame(model[!model[,1]=="#",],stringsAsFactors=FALSE)
}
