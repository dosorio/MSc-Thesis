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

# Lectura FPKM Astrocitos Saludables
Astrocyte_Expression <- read.csv("Data/GSE73721.csv",row.names = "Gene")[,c(9:26)]
Astrocyte_Activated <- Astrocyte_Expression[(rowSums(scale(Astrocyte_Expression)>=0)/ncol(Astrocyte_Expression))>=0.5,]

# Extraigo ID's genes en diferentes bases de datos
Human <- (UniProt.ws(taxId=9606))
Astrocyte_Genes <- select(Human,keys=toupper(rownames(Astrocyte_Activated)),columns = c("ENTREZ_GENE","ENSEMBL","EC"),keytype ="GENECARDS")

# Extrayendo de RECON
RECON <- read.csv("Data/RECON.csv",sep = "\t",stringsAsFactors = FALSE)
RECON$GPR <- gsub("([[:alnum:]]+)\\.[[:digit:]]+"," \\1 ",RECON$GPR)
RECON$GPR <- sapply(RECON$GPR, function(gpr){
  woSpaces <- gsub("\\(|\\)|[[:blank:]]+","",gpr)
  paste0(unique(unlist(lapply(lapply(unlist(strsplit(woSpaces,"or")), function(gpr){strsplit(gpr,"and")}),function(gpr){paste0("( ",paste0(sort(unlist(gpr)),collapse = " and ")," )")}))),collapse = " or ")
},USE.NAMES = FALSE)

# Extrayendo las reacciones asociadas a los genes en RECON
reactions <- sapply(unique(Astrocyte_Genes$ENTREZ_GENE[!is.na(Astrocyte_Genes$EC)]),function(enzyme){RECON$REACTION[grep(paste0("[[:blank:]]",enzyme,"[[:blank:]]"),RECON$GPR)]})
reactions <- unique(unlist(reactions))

# GapFind y GapFill
reactions <- gapFill(reactions,RECON$REACTION[RECON$GPR==""],consensus = TRUE)

# Añadiendo flujo
convert2sbml(RECON,"Results/DMEM.xml")
DMEM <- readSBMLmod("Results/DMEM.xml")

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
lowbnd(DMEM)[react_id(DMEM) == 'EX_pyr(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_thymd(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_tyr_L(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_o2(e)'] <- -0.530
lowbnd(DMEM)[react_id(DMEM) == 'EX_h2o(e)'] <- -100
lowbnd(DMEM)[react_id(DMEM) == 'EX_cl(e)'] <- -1000
lowbnd(DMEM)[react_id(DMEM) == 'EX_co2(e)'] <- 0.515
lowbnd(DMEM)[react_id(DMEM) == 'EX_so4(e)'] <- -100
lowbnd(DMEM)[react_id(DMEM) == 'EX_hdca(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_estradiol(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_tibolone(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_nh4(e)'] <- -100
uppbnd(DMEM)[react_id(DMEM)%in%react_id(findExchReact(DMEM))] <- 0
uppbnd(DMEM)[react_id(DMEM) == 'EX_prostgd2(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_prostge1(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_prostge2(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_prostgf2(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_prostgh2(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_prostgi2(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_HC02202(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_o2(e)'] <- -0.515
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
uppbnd(DMEM)[react_id(DMEM) == 'EX_co2(e)'] <- 0.530
uppbnd(DMEM)[react_id(DMEM) == 'EX_lac_L(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_glc_D(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_gln_L(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_ser_D(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_atp(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_taur(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_no(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_o2s(e)'] <- 1000
uppbnd(DMEM)[react_id(DMEM) == 'EX_h2o2(e)'] <- 1000

Reactions_Flux <- (RECON[getFluxDist(optimizeProb(DMEM))!=0,3])

DMEM@obj_coef <- rep(0,DMEM@react_num)

DMEM <- addReact(DMEM, id="MC", met=c("glu_L[e]","gln_L[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("gly[c]","ser_D[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("ser_L[c]","ser_D[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("glc_D[e]","lac_L[e]"),
                 Scoef=c(-1,2), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("glc_D[c]","atp[e]"),
                 Scoef=c(-1,2), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("cys_L[e]","gthrd[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","prostgd2[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","prostge1[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","prostge2[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","prostgf2[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","prostgh2[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","prostgi2[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","HC02202[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","HC02203[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","HC02204[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","HC02205[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","HC02206[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","HC02207[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","HC02208[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","HC02210[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","HC02213[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","HC02214[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","HC02216[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","HC02217[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","leuktrA4[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","leuktrB4[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","leuktrC4[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","leuktrD4[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","leuktrE4[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","leuktrF4[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","o2s[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","h2o2[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("arachd[c]","no[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

# 
DMEM <- rmReact(model = DMEM, react = "MC")

# Construyendo la reconstrucción
RECON$LOWER.BOUND <- DMEM@lowbnd
RECON$UPPER.BOUND <- DMEM@uppbnd
Astrocyte_Draft<- mapReactions(reactionList = unique(c(reactions,Reactions_Flux)),referenceData = RECON,by = "REACTION")
convert2sbml(Astrocyte_Draft,"Results/Astrocyte_Draft.xml")

#
Astrocyte_DraftM <- readSBMLmod("Results/Astrocyte_Draft.xml")
fv <- fluxVar(Astrocyte_DraftM)
t1 <- Astrocyte_DraftM@react_id[(fv@lp_obj[1:Astrocyte_DraftM@react_num]==0 & fv@lp_obj[(Astrocyte_DraftM@react_num+1):(2*Astrocyte_DraftM@react_num)]==0)]
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