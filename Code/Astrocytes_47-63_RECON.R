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

# Extraigo rutas metabolicas involucradas en la diferenciación de astrocito maduro
AstrocyteMetabolism<-enrichPathway(unique(na.omit(HSA$ENTREZ_GENE)),organism = "human",readable = TRUE,pvalueCutoff = 0.05)
ResumenAstrocyte<-summary(AstrocyteMetabolism)

# Extrayendo de HMR y RECON
RECON<-read.csv("Data/RECON_rxn.txt",sep = "\t")
RECON$Gene.reaction.association<-gsub("\\(","\\( ",RECON$Gene.reaction.association)
RECON$Gene.reaction.association<-gsub("\\)"," \\)",RECON$Gene.reaction.association)
#HMR <- read.csv("Data/rxn_HMRdatabase2_00.csv",sep = ";")
RXN<-NULL
for(i in HSA$ENTREZ_GENE){
  RXN<-unique(c(RXN,as.vector(RECON$Rxn.name[grep(paste0("[[:blank:]]",i,"\\.","[[:digit:]]+[[:blank:]]"),RECON$Gene.reaction.association)])))
}

# GAP2F
as.searchable <- function(metabolite){
  metabolite <- gsub("\\(","\\\\(",metabolite)
  metabolite <- gsub("\\)","\\\\)",metabolite)
  metabolite <- gsub("\\[","\\\\[",metabolite)
  metabolite <- gsub("\\]","\\\\]",metabolite)
  metabolite <- gsub("\\-","\\\\-",metabolite)
  metabolite <- gsub("\\+","\\\\+",metabolite)
  metabolite <- gsub("\\_","\\\\_",metabolite)
  return(metabolite)
}

reactant_fill <- function(metabolite,reference){
  unique(c(as.vector(reference)[grep(paste0("<?=>?[[:print:]]+?",as.searchable(metabolite)),reference)],
           as.vector(reference)[grep(paste0("<=>[[:print:]]+?",as.searchable(metabolite)),reference)],
           as.vector(reference)[grep(paste0(as.searchable(metabolite),"[[:print:]]+?<=>"),reference)]))
}
product_fill <- function(metabolite,reference){
  unique(c(as.vector(reference)[grep(paste0(as.searchable(metabolite),"[[:print:]]+?<?=>?"),reference)],
           as.vector(reference)[grep(paste0("<=>[[:print:]]+?",as.searchable(metabolite)),reference)],
           as.vector(reference)[grep(paste0(as.searchable(metabolite),"[[:print:]]+?<=>"),reference)]))
}


fill <- NULL
reactionList<-unique(c(as.vector(RECON$Formula[RECON$Rxn.name%in%RXN]),fill))
repeat{
  ometabolites <- (orphan.reactants(reactionList))
  #ometabolites<- ometabolites[-grep("\\[x\\]",ometabolites)]
  fill <- as.vector(unlist(sapply(ometabolites, function(x){reactant_fill(x,RECON$Formula[!RECON$Formula%in%reactionList])},simplify = TRUE)))
  reactionList<-unique(as.vector(c(reactionList,fill)))
  if (length(fill)==0){
    break
  }
}
# Reconstrucción Borrador
Astrocyte_DRAFT<-RECON[RECON$Formula%in%reactionList,]

orphan<- unique(orphan.reactants(as.vector(Astrocyte_DRAFT$Formula)))


# Convirtiendo a Sybil
cobra2sybil<-function(reaction){
  coeficients <- function(met) {
    coef<-regmatches(met, gregexpr('^[[:digit:]][[:graph:]]*[[:blank:]]', met))
    gsub("[[:blank:]]","",coef)
  }
  
  if(grepl("[[:blank:]]=>[[:blank:]]",reaction)){
    reaction <-unlist(strsplit(reaction,"=>"))
    reactant <- reaction[1]
    reactant <- unlist(strsplit(reactant,"[[:blank:]]\\+[[:blank:]]"))
    reactant <- gsub("^[[:blank:]]*","",reactant)
    reactant <- gsub("[[:blank:]]*$","",reactant)
    products <- reaction[2]
    products <- unlist(strsplit(products,"[[:blank:]]\\+[[:blank:]]"))
    products <- gsub("^[[:blank:]]*","",products)
    products <- gsub("[[:blank:]]*$","",products)
    r_coefic <- as.numeric(sapply(reactant, coeficients))
    r_coefic[is.na(r_coefic)]<-1
    p_coefic <- as.numeric(sapply(products, coeficients))
    p_coefic[is.na(p_coefic)]<-1
    reactant <- gsub('^[[:digit:]]([[:print:]]*)[[:blank:]]',"",reactant)
    products <- gsub('^[[:digit:]][[:graph:]]*[[:blank:]]',"",products)
    reactant <- mapply(function(c,m){paste("(",c,") ", m, sep = "")}, c=r_coefic, m=reactant)
    products <- mapply(function(c,m){paste("(",c,") ", m, sep = "")}, c=p_coefic, m=products)
    reaction <- paste(paste0(reactant,collapse = " + "), paste0(products,collapse = " + "),sep = " --> ")
  }else{
    reaction <-unlist(strsplit(reaction,"<=>"))
    reactant <- reaction[1]
    reactant <- unlist(strsplit(reactant,"[[:blank:]]\\+[[:blank:]]"))
    reactant <- gsub("^[[:blank:]]*","",reactant)
    reactant <- gsub("[[:blank:]]*$","",reactant)
    products <- reaction[2]
    products <- unlist(strsplit(products,"[[:blank:]]\\+[[:blank:]]"))
    products <- gsub("^[[:blank:]]*","",products)
    products <- gsub("[[:blank:]]*$","",products)
    r_coefic <- as.numeric(sapply(reactant, coeficients))
    r_coefic[is.na(r_coefic)]<-1
    p_coefic <- as.numeric(sapply(products, coeficients))
    p_coefic[is.na(p_coefic)]<-1
    reactant <- gsub('^[[:digit:]]([[:print:]]*)[[:blank:]]',"",reactant)
    products <- gsub('^[[:digit:]][[:graph:]]*[[:blank:]]',"",products)
    reactant <- mapply(function(c,m){paste("(",c,") ", m, sep = "")}, c=r_coefic, m=reactant)
    products <- mapply(function(c,m){paste("(",c,") ", m, sep = "")}, c=p_coefic, m=products)
    reaction <- paste(paste0(reactant,collapse = " + "), paste0(products,collapse = " + "),sep = " <==> ")
  }
  return(reaction)
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
  reaction <- gsub("\\([[:graph:]]*\\)[[:blank:]]+","",reaction)
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
orphans <- orphan.reactants(as.vector(Astrocyte_DRAFT$Formula))[!grepl("\\[x\\]",orphan.reactants(as.vector(Astrocyte_DRAFT$Formula)))]
reactions<-sapply(c(as.vector(Astrocyte_DRAFT$Formula)),cobra2sybil)
reactions<-as.vector(gsub("\\[x\\]","\\[b\\]",reactions))
names <- Astrocyte_DRAFT$Rxn.name
rever <- as.vector(sapply(reactions,reversible))
comp<- as.vector(sapply(reactions,compartment))
bound<-(!comp=="b")

Astrocyte_MODEL<-as.data.frame(cbind(as.character(names[bound]),as.character(names[bound]),as.character(reactions[bound]),as.character(rever[bound]),as.character(comp[bound]),-1000,1000,0,as.character(Astrocyte_DRAFT$Gene.reaction.association[bound]),""),stringsAsFactors = FALSE)

for (orphan in orphans){
  Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),] <- (cbind(paste0(gsub("\\[[[:alpha:]]\\]","[b]",orphan)),paste0("Ex_",orphan),paste0(gsub("\\[[[:alpha:]]\\]","[b]",orphan)," --> ",orphan),"irreversible","b",-1000,1000,0,"","Exchange"))
}



Astrocyte_MODEL$V6[rever[bound]=="irreversible"]<-0
#Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c("Biomass","OBJECTIVE1","glc_D[e] --> adp[c]","irreversible","c",0,1000,1,"","")
#Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c("Biomass","OBJECTIVE1","(20.7045) atp[c] + (20.6508) h2o[c] + (0.38587) glu_L[c] + (0.35261) asp_L[c] + (0.036117) gtp[c] + (0.27942) asn_L[c] + (0.50563) ala_L[c] + (0.046571) cys_L[c] + (0.326) gln_L[c] + (0.53889) gly[c] + (0.39253) ser_L[c] + (0.31269) thr_L[c] + (0.59211) lys_L[c] + (0.35926) arg_L[c] + (0.15302) met_L[c] + (0.023315) pail_hs[c] + (0.039036) ctp[c] + (0.15446) pchol_hs[c] + (0.055374) pe_hs[c] + (0.020401) chsterol[c] + (0.002914) pglyc_hs[c] + (0.011658) clpn_hs[c] + (0.053446) utp[c] + (0.009898) dgtp[n] + (0.009442) dctp[n] + (0.013183) datp[n] + (0.013091) dttp[n] + (0.27519) g6p[c] + (0.12641) his_L[c] + (0.15967) tyr_L[c] + (0.28608) ile_L[c] + (0.54554) leu_L[c] + (0.013306) trp_L[c] + (0.25947) phe_L[c] + (0.41248) pro_L[c] + (0.005829) ps_hs[c] + (0.017486) sphmyln_hs[c] + (0.35261) val_L[c] --> (20.6508) adp[c] + (20.6508) h[c] + (20.6508) pi[c]","irreversible","c",0,1000,1,"","")
Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c("Biomass","OBJECTIVE1","(20.7045) atp[c]  --> (20.6508) adp[c] + (20.6508) h[c] + (20.6508) pi[c]","irreversible","c",0,1000,1,"","")
#Astrocyte_MODEL[(dim(Astrocyte_MODEL)[1]+1),]<-c("GBM","OBJECTIVE2","--> oaa[m] + succ[m] + gthrd[c] + r5p[c]","irreversible","c",0,1000,1,"","")

colnames(Astrocyte_MODEL)<-c("abbreviation","name","equation","reversible","compartment","lowbnd","uppbnd","obj_coef","rule","subsystem")
write.table(Astrocyte_MODEL,sep = "\t",row.names = FALSE,file = "Astrocyte_react.tsv")

metabolites <- unique(unlist(sapply(Astrocyte_MODEL$equation, metabolites_f)))
m_names <- unique(gsub('\\[[[:alpha:]]+\\]$',"",metabolites))
compart <- as.vector(sapply(m_names, metmodel))
Astrocyte_METS <- as.data.frame(cbind(m_names,m_names,compart))
Astrocyte_METS <- Astrocyte_METS[!Astrocyte_METS$abbreviation=="",]
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


