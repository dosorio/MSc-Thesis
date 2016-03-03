setwd("~/Dropbox/Maestría Bioinformática/Tesis de Maestría/")
source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq")
#biocLite("ReactomePA")
#biocLite("KEGGREST")
#biocLite("reactome.db")
#biocLite("UniProt.ws")
#biocLite("org.Hs.eg.db")
#install_github("dosorio/convert2sybil")
# require(DESeq)
# require(edgeR)
require(caret)
# require(KEGGREST)
require(reactome.db)
require(UniProt.ws)
require(ReactomePA)

# Lectura FPKM astrocitos 21 - 63 yo
Astrocyte <-read.csv("Data/GSE73721.csv",row.names = "Gene")[,c(18:26)]

# Eliminando genes que no varian expresión en ningun tiempo (ruido basal)
ActiveGenes<-Astrocyte[-nearZeroVar(t(Astrocyte)),]

# Basandome en el FPKM (anotando nuevamente), busco cuales son diferencialmente expresados (mayor expresión en >40) usando prueba t
DE<-NULL
for(i in 1:dim(ActiveGenes)[1]){
  DE[i]<-try(t.test(t(ActiveGenes)[,i]~c(rep('A',3),rep('B',6)), alternative="less")$p.value,silent = TRUE)
}

# Identifico genes diferencialmente expresados usando p<0.5
DE<-as.numeric(DE)<0.05
table(DE)

# Extraigo ID's genes en diferentes bases de datos
Human<-(UniProt.ws(taxId=9606))
HSA<-select(Human,keys=toupper(unique(rownames(ActiveGenes))),columns = c("ENTREZ_GENE","EC"),keytype ="GENECARDS")
#MatureAstrocyte<- select(Human,keys=toupper(unique(rownames(ActiveGenes[DE,]))),columns = c("ENSEMBL","ENTREZ_GENE","ENTRY-NAME","KO","KEGG","EC","REACTOME","PROTEIN-NAMES"),keytype ="GENECARDS")

# Excluyo proteínas con actividad enzimática
EC<-HSA[!is.na(HSA$EC),]

# Extraigo rutas metabolicas involucradas en la diferenciación de astrocito maduro
AstrocyteMetabolism<-enrichPathway(unique(na.omit(HSA$ENTREZ_GENE)),organism = "human",readable = TRUE,pvalueCutoff = 0.05)
ResumenAstrocyte<-summary(AstrocyteMetabolism)

# Extrayendo ENSEMBL ID's
EC_C<-unique(unlist(strsplit(EC$EC,";")))

# Extrayendo de HMR
setwd("~/Dropbox/Maestría Bioinformática/Tesis de Maestría/")
RECON<-read.csv("Data/RECON_rxn.txt",sep = "\t")
RXN<-NULL
for (i in 1:length(EC$ENTREZ_GENE)){
  RXN<-unique(c(RXN,as.vector(RECON[grep(paste0("[[:punct:]]",EC$ENTREZ_GENE[i],"\\."),RECON$Gene.reaction.association),1])))
}

# Reconstrucción Borrador
Astrocyte_DRAFT<-RECON[RECON$Rxn.name%in%RXN,]

# Convirtiendo a Sybil
cobra2sybil<-function(reaction){
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

reversible <- function(reaction){
  if(grepl("[[:blank:]]-->[[:blank:]]",reaction)){
    return("irreversible")
  } else {
    return("reversible")
  }
}

compartment <- function(reaction){
  if(grepl("[[:blank:]]=>[[:blank:]]",reaction)){
    reaction <- strsplit(reaction," <==> ",fixed = TRUE)[[1]]
  } else {
    reaction <- strsplit(reaction," --> ",fixed = TRUE)[[1]]}
  reaction <- unlist(strsplit(reaction," + ",fixed = TRUE))
  comp<-unlist(regmatches(reaction, gregexpr('\\[[[:alpha:]]\\]$', reaction)))
  comp<-gsub("\\[","",comp)
  comp<-gsub("\\]","",comp)
  paste(unique(comp),collapse = ", ")
}

metabolites_f<-function(reaction){
  if(grepl("[[:blank:]]?-->[[:blank:]]?",reaction)){
    reaction <- strsplit(reaction,"[[:blank:]]?-->[[:blank:]]?")[[1]]
  } else {
    reaction <- strsplit(reaction,"[[:blank:]]?<==>[[:blank:]]?")[[1]]}
  reaction <- unlist(strsplit(reaction,"[[:blank:]]?\\+[[:blank:]]?"))
  reaction <- gsub("\\([[:digit:]]+\\)[[:blank:]]","",reaction)
  gsub("[[:blank:]]","",reaction)
}

metmodel<-function(m_name){
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

setwd("~/Desktop/")
reactions<-sapply(Astrocyte_DRAFT$Formula,cobra2sybil)
names <- Astrocyte_DRAFT$Rxn.name
rever <- as.vector(sapply(reactions,reversible))
comp<- as.vector(sapply(reactions,compartment))

Astrocyte_MODEL<-as.data.frame(cbind(as.character(names),as.character(Astrocyte_DRAFT$Rxn.description),as.character(reactions),as.character(rever),as.character(comp),-1000,1000,0,as.character(Astrocyte_DRAFT$Gene.reaction.association),""),stringsAsFactors = FALSE)
Astrocyte_MODEL$V6[rever=="irreversible"]<-0
m_out <- c("L-lactate[e]",
           "glutamine[e]",
           "sphingosine-1-phosphate[e]",
           "starch structure 1[c]",
            "glycogenin G11[c]",
            "prostaglandin E2[e]",
            "prostaglandin F2alpha[e]",
            "prostaglandin E1[e]",
            "5,6-EET[c]",
            "8,9-EET[c]",
            "11,12-EET[c]",
            "14,15-EET[c]",
            "GSH[e]",
           "ATP[m]",
           "glycogen[c]",
           "H2O[c]",
           "H[c]",
           "CO2[c]",
           "glucose[e]",
           "OAA[m]",
           "succinate[m]",
           "ribose-5-phosphate[c]"
)
m_in <- c("glucose[e]",
          "glutamate[e]",
          "O2[c]",
          "Pi[e]",
          "NH3[e]",
          "H[e]",
          "Na[e]",
          "Fe2[e]",
          "Li[e]",
          "Ca2[e]",
          "Mg2[e]",
          "Cu2[e]",
          "K[e]",
          "folate[e]",
          "NO[e]",
          "NH4[e]",
          "zinc[e]"
          # "carbonate[e]",
          # "cysteine[e]",
          # "glycine[e]",
          # "arginine[e]",
          # "histidine[e]",
          # "methionine[e]",
          # "tyrosine[e]",
          # "glutamine[e]",
          # "serine[e]",
          # "tryptophan[e]",
          # "phenylalanine[e]",
          # "leucine[e]",
          # "valine[e]",
          # "threonine[e]",
          # "isoleucine[e]",
          # "cystine[e]"
)

for(compound in m_out){
  id = paste0(compound,"OUT")
  Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c(id,id,paste0(compound," --> "),"irreversible","e",0,1000,0,"","")
}
for(compound in m_in){
  id = paste0(compound,"IN")
  Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c(id,id,paste0(compound,"--> "),"irreversible","e",-1000,0,0,"","")
}

Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c("GLIOT","OBJECTIVE1","beta-D-glucose[c] + glutamate[e] + O2[c] + Pi[m] + H[m] --> H2O[m] + ATP[m] + H[m] + OAA[m] + succinate[m] + ribose-5-phosphate[c] + glycogenin G11[c] + L-lactate[e] + glutamine[e] + sphingosine-1-phosphate[e] + prostaglandin E2[e] + prostaglandin F2alpha[e] + prostaglandin E1[e] + 5,6-EET[c] + 8,9-EET[c] + 11,12-EET[c] + 14,15-EET[c] + GSH[e]","irreversible","c",0,1000,1,"","")
#Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c("ATPSYN","OBJECTIVE2","ADP[m] + Pi[m] + H[m] --> H2O[m] + ATP[m] + H[m]","irreversible","c",0,1000,1,"","")
#Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c("GBM","OBJECTIVE3","--> OAA[m] + succinate[m] + GSH[c] + ribose-5-phosphate[c]","irreversible","c",0,1000,1,"","")

colnames(Astrocyte_MODEL)<-c("abbreviation","name","equation","reversible","compartment","lowbnd","uppbnd","obj_coef","rule","subsystem")
write.table(Astrocyte_MODEL,sep = "\t",row.names = FALSE,file = "Astrocyte_react.tsv")


<<<<<<< HEAD
metabolites <- unique(unlist(sapply(as.vector(Astrocyte_MODEL$equation),metabolites_f)))
m_names <- unique(gsub('\\[[[:alpha:]]+\\]$',"",metabolites))
compart <- as.vector(sapply(m_names, metmodel))
Astrocyte_METS <- as.data.frame(cbind(m_names,m_names,compart))
colnames(Astrocyte_METS)<-c("abbreviation","name","compartment")
write.table(Astrocyte_METS,sep = "\t",row.names = FALSE,file = "Astrocyte_met.tsv")

name = "T1"
id = "T1"
description = ""
abbreviation = ""
Nmetabolites = ""
Nreactions = ""
Ngenes = ""
Nnnz = ""
Astrocyte_DESC<-cbind(name,id)
write.table(Astrocyte_DESC,sep = "\t",row.names = FALSE,file = "Astrocyte_desc.tsv")
