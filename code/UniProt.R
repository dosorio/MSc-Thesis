# UniProt.R
# Retrieve EC numbers from UniProt
# Daniel Camilo Osorio
# Maestría en Bioinformática - Universidad Nacional de Colombia (Bogotá)
# Laboratorio de Bioquímica Teórica y Bioinformática - Pontificia Universidad Javeriana (Bogotá)

# Setting the working directory
setwd("Documents/Proteomics/")

# Loading the UniProt.ws package
library(UniProt.ws)

# Reading the file
Proteome<-read.table(commandArgs(TRUE)[1],sep = "\t")
#Proteome<-read.table("results/Human_OligoDendocyte/OL-XO_Proteins.txt",sep = "\t")

# Filtering the confident proteins
Proteome<-Proteome[which(Proteome$V22 == "Confident"),]

# Replacing the protein list
Specie<-(UniProt.ws(taxId=9606))
Specie@taxIdUniprots<- as.vector(Proteome$V2)

# Retrieving the UNIPROT EC numbers
Data <- select(Specie, keys(Specie,"UNIPROTKB"), "EC", "UNIPROTKB")

# Writing output files
write.table(Data[!is.na(Data$EC),],"results/Human_OligoDendocyte/EC_OL-OX.txt",sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
write.table(Data[is.na(Data$EC),],"results/Human_OligoDendocyte/NoEC_OL-OX.txt",sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
