setwd("~/Dropbox/Maestria Bioinformatica/Tesis de Maestria/")
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
require(minval)

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
HSA<-select(Human,keys=toupper(unique(rownames(ActiveGenes))),columns = c("ENTREZ_GENE","ENSEMBL","EC"),keytype ="GENECARDS")
#MatureAstrocyte<- select(Human,keys=toupper(unique(rownames(ActiveGenes[DE,]))),columns = c("ENSEMBL","ENTREZ_GENE","ENTRY-NAME","KO","KEGG","EC","REACTOME","PROTEIN-NAMES"),keytype ="GENECARDS")

# Excluyo proteínas con actividad enzimática
EC<-HSA[!is.na(HSA$EC),]

# Extraigo rutas metabolicas involucradas en la diferenciación de astrocito maduro
AstrocyteMetabolism<-enrichPathway(unique(na.omit(HSA$ENTREZ_GENE)),organism = "human",readable = TRUE,pvalueCutoff = 0.05)
ResumenAstrocyte<-summary(AstrocyteMetabolism)

# Extrayendo ENSEMBL ID's
EC_C<-unique(unlist(strsplit(EC$EC,";")))

# Extrayendo de HMR
HMR<-read.csv("Data/rxn_HMRdatabase2_00.csv",sep = ";")
RXN<-NULL
for (i in 1:length(EC$ENTREZ_GENE)){
  RXN<-unique(c(RXN,as.vector(HMR[grep(EC$ENSEMBL[i],HMR$GENE.ASSOCIATION),2])))
}

# Reconstrucción Borrador
Astrocyte_DRAFT<-HMR[HMR$RXNID%in%RXN,]
orphan<-unique(orphan.reactants(Astrocyte_DRAFT$EQUATION))



# Convirtiendo a Sybil
cobra2sybil<-function(reaction){
  coeficients <- function(met) {
    coef<-regmatches(met, gregexpr('^[[:digit:]]+[[:punct:]]?[[:digit:]]?[[:digit:]]?[[:blank:]]', met))
    gsub("[[:blank:]]","",coef)
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
    reactant <- gsub('^[[:digit:]]+[[:punct:]]?[[:digit:]]?[[:digit:]]?[[:blank:]]',"",reactant)
    products <- gsub('^[[:digit:]]+[[:punct:]]?[[:digit:]]?[[:digit:]]?[[:blank:]]',"",products)
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
  if(grepl("[[:blank:]]<==>[[:blank:]]",reaction)){
    reaction <- unlist(strsplit(reaction," <==> ",fixed = TRUE))
  } else {
    reaction <- unlist(strsplit(reaction," --> ",fixed = TRUE))}
  reaction <- unlist(strsplit(reaction," + ",fixed = TRUE))
  comp<-unlist(regmatches(reaction, gregexpr('\\[[[:alpha:]]\\]$', reaction)))
  comp<-gsub("\\[","",comp)
  comp<-gsub("\\]","",comp)
  paste(unique(comp),collapse = ", ")
}

metabolites_f<-function(reaction){
  if(grepl("[[:blank:]]<==>[[:blank:]]",reaction)){
    reaction <- unlist(strsplit(reaction," <==> ",fixed = TRUE))
  } else {
    reaction <- unlist(strsplit(reaction," --> ",fixed = TRUE))}
  reaction <- unlist(strsplit(reaction," + ",fixed = TRUE))
  reaction <- gsub("\\([[:digit:]]+[[:punct:]]?[[:digit:]]?[[:digit:]]?\\)[[:blank:]]","",reaction)
  reaction
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
reactions<-sapply(Astrocyte_DRAFT$EQUATION,cobra2sybil)
names <- Astrocyte_DRAFT$RXNID
rever <- as.vector(sapply(reactions,reversible))
comp<- as.vector(sapply(reactions,compartment))

Astrocyte_MODEL<-as.data.frame(cbind(as.character(names),as.character(names),as.character(reactions),as.character(rever),as.character(comp),-1000,1000,0,"",""),stringsAsFactors = FALSE)
Astrocyte_MODEL$V6[rever=="irreversible"]<-0
m_out <- c("ATP[m]","H[m]")
m_in <- unique(c(gsub("\\[s\\]","\\[e\\]",gsub("\\+","",orphan)),"ADP[m]","Pi[m]","H[m]"))

for(compound in m_out){
  id = paste0(compound,"OUT")
  Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c(id,id,paste0("(1) ",compound," --> "),"irreversible","e",0,1000,0,"","")
}
for(compound in m_in){
  id = paste0(compound,"IN")
  Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c(id,id,paste0("(1) ",compound," --> "),"irreversible","e",-1000,0,0,"","")
}

Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c("GLIOT","OBJECTIVE1","ADP[m] + Pi[m] + H[m] --> ATP[m] + H[m]","irreversible","c",0,1000,1,"","")
#Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c("ATPSYN","OBJECTIVE2","ADP[m] + Pi[m] + H[m] --> H2O[m] + ATP[m] + H[m]","irreversible","c",0,1000,1,"","")
#Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c("GBM","OBJECTIVE3","--> OAA[m] + succinate[m] + GSH[c] + ribose-5-phosphate[c]","irreversible","c",0,1000,1,"","")

colnames(Astrocyte_MODEL)<-c("abbreviation","name","equation","reversible","compartment","lowbnd","uppbnd","obj_coef","rule","subsystem")
write.table(Astrocyte_MODEL,sep = "\t",row.names = FALSE,file = "Astrocyte_react.tsv")

metabolites <- unique(unlist(sapply(Astrocyte_MODEL$equation, metabolites_f)))
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

