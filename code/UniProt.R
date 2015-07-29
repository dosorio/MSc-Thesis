# UniProt.R
# Retrieve EC numbers from UniProt
# Daniel Camilo Osorio
# Maestría en Bioinformática - Universidad Nacional de Colombia (Bogotá)
# Laboratorio de Bioquímica Teórica y Bioinformática - Pontificia Universidad Javeriana (Bogotá)

# Load the UniProt.ws package
library(UniProt.ws)

# Read the file
Proteome<-read.table(commandArgs(TRUE)[1],sep = "\t")
#Proteome<-read.table("Documents/Proteomics/results/Human_OligoDendocyte/OL-XO_Proteins.txt",sep = "\t")

# Filter the confident proteins
Proteome<-Proteome[which(Proteome$V22 == "Confident"),]

# Replace the protein list
Specie<-(UniProt.ws(taxId=9606))
Specie@taxIdUniprots<- as.vector(Proteome$V2)


ID<-keys(Specie,"UNIPROTKB")
keys <- ID
columns <- c("EC")
kt <- c("UNIPROTKB")
Data <- select(Specie, keys(Specie,"UNIPROTKB"), "EC", "UNIPROTKB")
NoEC<-Data[is.na(Data$EC),]
EC<-Data[!is.na(Data$EC),]
write.table(EC,"EC_OL-OX.txt",sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
write.table(NoEC,"NoEC_OL-OX.txt",sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)


# RXN<-NULL
# for (EC in paste("EC:",IDs,sep = "")){
#   RXN<-append(RXN,(which(HMA$V5 == EC )))
# }

