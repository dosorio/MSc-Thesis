library(UniProt.ws)
Proteome<-read.table("Documents/Proteomics/Results/Human_OligoDendocyte/OL_Proteins.txt",sep = "\t")
Proteome<-Proteome[which(Proteome$V22 == "Confident"),]
Specie<-(UniProt.ws(taxId=9606))
Specie@taxIdUniprots<- as.vector(Proteome$V2)

#columns(Specie)


ID<-keys(Specie,"UNIPROTKB")
keys <- ID
columns <- c("EC")
kt <- c("UNIPROTKB")
Data <- select(Specie, keys, columns, kt)
EC<-Data[!is.na(Data$EC),]

