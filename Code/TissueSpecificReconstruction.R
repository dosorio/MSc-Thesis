# install.packages("devtools")
# source("https://bioconductor.org/biocLite.R")
# library(devtools)
# biocLite("UniProt.ws")
# biocLite("GEOquery")
# install_github("gibbslab/g2f")
# install_github("gibbslab/minval")
# biocLite("ReactomePA")
# install.packages("sybilSBML")
require(UniProt.ws)
require(GEOquery)
require(g2f)
require(minval)
require(ReactomePA)
require(sybilSBML)

# Descarga FPKM Astrocitos de Corteza Saludables
getGEOSuppFiles(GEO = "GSE73721",makeDirectory = TRUE,baseDir = "Data/")

# Removiendo otros tipos de datos
Astrocyte_Expression <- read.csv(file = "~/Documents/masterThesis/Data/GSE73721/GSE73721_Human_and_mouse_table.csv.gz",
                                 header = TRUE,
                                 row.names = "Gene")[,c(9:26)]

# Identificando genes expresados por encima de la media en al menos el 50% de los casos
Astrocyte_Activated <- Astrocyte_Expression[(rowSums(scale(Astrocyte_Expression)>=0)/ncol(Astrocyte_Expression))>=0.5,]

# Extraigo ID's genes en diferentes bases de datos
Human <- (UniProt.ws(taxId=9606))
Astrocyte_Genes <- select(Human,keys=toupper(rownames(Astrocyte_Activated)),columns = c("ENTREZ_GENE","ENSEMBL","EC"),keytype ="GENECARDS")

# Leyendo RECON para usarla como referencia
RECON <- read.csv(file = "Data/RECON.csv",
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  header = TRUE,
                  col.names = c("ID","DESCRIPTION","REACTION","GPR","REVERSIBLE","LOWER.BOUND","UPPER.BOUND","OBJECTIVE"))

# Convirtiendo las GPR a ENTREZ
RECON$GPR <- gsub("([[:alnum:]]+)\\.[[:digit:]]+"," \\1 ",RECON$GPR)
RECON$REACTION <- gsub("[[:blank:]]+"," ",RECON$REACTION)
RECON$GPR <- sapply(RECON$GPR, function(gpr){
  woSpaces <- gsub("\\(|\\)|[[:blank:]]+","",gpr)
  paste0(unique(unlist(lapply(lapply(unlist(strsplit(woSpaces,"or")), function(gpr){strsplit(gpr,"and")}),function(gpr){paste0("( ",paste0(sort(unlist(gpr)),collapse = " and ")," )")}))),collapse = " or ")
},USE.NAMES = FALSE)

# Extrayendo las reacciones asociadas a los genes en RECON
Astrocyte_Reactions <- unique(unlist(sapply(unique(Astrocyte_Genes$ENTREZ_GENE[!is.na(Astrocyte_Genes$EC)]),function(enzyme){RECON$REACTION[grep(paste0("[[:blank:]]",enzyme,"[[:blank:]]"),RECON$GPR)]})))

# GapFind y GapFill
Astrocyte_Reactions <- gapFill(reactionList = Astrocyte_Reactions,
                               reference = RECON$REACTION[nchar(RECON$GPR)==0],
                               consensus = TRUE)

# Añadiendo flujo
convert2sbml(RECON,"Results/RECON.xml")
DMEM <- readSBMLmod("Results/RECON.xml")

# DMEM
lowbnd(DMEM)[react_id(DMEM)%in%react_id(findExchReact(DMEM))] <- 0
lowbnd(DMEM)[react_id(DMEM) == 'EX_ca2(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_glc(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_fe3(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_k(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_na1(e)'] <- -1000
lowbnd(DMEM)[react_id(DMEM) == 'EX_HC02172(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_Rtotal(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_Rtotal2(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_Rtotal3(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_pi(e)'] <- -100
lowbnd(DMEM)[react_id(DMEM) == 'EX_ala_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_arg_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_asn_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_asp_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_cys_L(e)'] <- -1000
lowbnd(DMEM)[react_id(DMEM) == 'EX_glu_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_gly(e)'] <- -1000
lowbnd(DMEM)[react_id(DMEM) == 'EX_his_L(e)'] <- -1000
lowbnd(DMEM)[react_id(DMEM) == 'EX_ile_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_leu_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_lys_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_met_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_phe_L(e)'] <- -1000
lowbnd(DMEM)[react_id(DMEM) == 'EX_pro_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_ser_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_thr_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_trp_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_val_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_btn(e)'] <- -100
lowbnd(DMEM)[react_id(DMEM) == 'EX_chol(e)'] <- -1000
lowbnd(DMEM)[react_id(DMEM) == 'EX_aqcobal(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_fol(e)'] <- -1000
lowbnd(DMEM)[react_id(DMEM) == 'EX_inost(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_lnlc(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_lnlnca(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_ncam(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_pydxn(e)'] <- -1000
lowbnd(DMEM)[react_id(DMEM) == 'EX_ribflv(e)'] <- -1000
lowbnd(DMEM)[react_id(DMEM) == 'EX_thm(e)'] <- -1000
lowbnd(DMEM)[react_id(DMEM) == 'EX_thymd(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_tyr_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_o2(e)'] <- -0.530
lowbnd(DMEM)[react_id(DMEM) == 'EX_h2o(e)'] <- -100
lowbnd(DMEM)[react_id(DMEM) == 'EX_nh4(e)'] <- -100
lowbnd(DMEM)[react_id(DMEM) == 'EX_cl(e)'] <- -1000
lowbnd(DMEM)[react_id(DMEM) == 'EX_co2(e)'] <- 0.515
lowbnd(DMEM)[react_id(DMEM) == 'EX_so4(e)'] <- -100
lowbnd(DMEM)[react_id(DMEM) == 'EX_hdca(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_estradiol(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_nh4(e)'] <- -100
lowbnd(DMEM)[react_id(DMEM) == 'EX_4abut(e)'] <- -1
uppbnd(DMEM)[react_id(DMEM)%in%react_id(findExchReact(DMEM))] <- 0
uppbnd(DMEM)[react_id(DMEM) == 'EX_co2(e)'] <- 0.530
uppbnd(DMEM)[react_id(DMEM) == 'EX_o2(e)'] <- -0.515
uppbnd(DMEM)[react_id(DMEM) == 'EX_prostgd2(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_prostge1(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_prostge2(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_prostgf2(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_prostgh2(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_prostgi2(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_HC02202(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_HC02203(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_HC02204(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_HC02205(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_HC02206(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_HC02207(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_HC02208(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_HC02210(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_HC02213(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_HC02214(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_HC02216(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_HC02217(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_leuktrA4(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_leuktrB4(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_leuktrC4(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_leuktrD4(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_leuktrE4(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_leuktrF4(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_gthrd(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_lac_L(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_glc_D(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_glu_L(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_gln_L(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_ser_D(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_atp(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_no(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_o2s(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_fe2(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_h2o2(e)'] <- 1000


Enrichment <- (RECON[getFluxDist(optimizeProb(DMEM))!=0,3])

DMEM@obj_coef <- rep(0,DMEM@react_num)

DMEM <- addReact(DMEM, id="MC", met=c("glu_L[e]","gln_L[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("nh4[c]","glu_L[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("gly[e]","ser_D[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("ser_L[c]","ser_D[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("glc_D[e]","lac_L[e]"),
                 Scoef=c(-1,2), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("glc_D[c]","atp[e]"),
                 Scoef=c(-1,2), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("cys_L[e]","gthrd[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","prostgd2[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","prostge1[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","prostge2[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","prostgh2[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","prostgi2[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","HC02208[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","HC02210[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","HC02213[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","HC02214[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","HC02216[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","HC02217[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","leuktrB4[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","leuktrE4[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("hdca[c]","no[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Enrichment <- unique(c(Enrichment,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

# 
DMEM <- rmReact(model = DMEM, react = "MC")

# Construyendo la reconstrucción
RECON$LOWER.BOUND <- DMEM@lowbnd
RECON$UPPER.BOUND <- DMEM@uppbnd
Astrocyte_Draft<- mapReactions(reactionList = unique(c(Astrocyte_Reactions,Enrichment)),
                               referenceData = RECON,
                               by = "REACTION")
convert2sbml(Astrocyte_Draft,"Results/Astrocyte_Draft.xml")

#
Astrocyte_DraftM <- readSBMLmod("Results/Astrocyte_Draft.xml")
optimizeProb(Astrocyte_DraftM)

# 
woFlux <- blockedReactions(Astrocyte_DraftM)

# 
Astrocyte_Reconstruction <- mapReactions(reactionList = woFlux,
                                         referenceData = Astrocyte_Draft,
                                         by = "ID",
                                         inverse = TRUE)

#
convert2sbml(Astrocyte_Reconstruction,"Results/Astrocyte.xml")

#
Astrocyte <- readSBMLmod("Results/Astrocyte.xml")
optimizeProb(Astrocyte)

# FBA
model_FBA <- optimizeProb(Astrocyte)
model_FBA

# Minimizing the total absolute fluxes MTF (Evaluando las múltiples posibles soluciones)
model_MTF <- optimizeProb(Astrocyte, algorithm = "mtf", wtobj = mod_obj(model_FBA))
getNetFlux(getFluxDist(model_MTF,findExchReact(Astrocyte)))
