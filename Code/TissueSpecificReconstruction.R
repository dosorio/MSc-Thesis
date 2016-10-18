# install.packages("devtools")
# source("https://bioconductor.org/biocLite.R")
# library(devtools)
# biocLite("UniProt.ws")
# install_github("gibbslab/g2f")
# install_github("gibbslab/minval")
# biocLite("ReactomePA")
# install.packages("sybilSBML")
require(UniProt.ws)
require(g2f)
require(minval)
require(ReactomePA)
require(sybilSBML)

# Lectura FPKM Astrocitos Saludables
Astrocyte_Expression <- read.csv("Data/GSE73721.csv",row.names = "Gene")[,c(9:26)]
Astrocyte_Activated <- Astrocyte_Expression[(rowSums(scale(Astrocyte_Expression)>0)/ncol(Astrocyte_Expression))>=0.5,]

# Extraigo ID's genes en diferentes bases de datos
Human <- (UniProt.ws(taxId=9606))
Astrocyte_Genes <- select(Human,keys=toupper(rownames(Astrocyte_Activated)),columns = c("ENTREZ_GENE","ENSEMBL","EC"),keytype ="GENECARDS")

# Extrayendo de RECON
RECON <- as.data.frame.array(read.csv("Data/RECON_rxn.txt",sep = "\t"))
RECON$Gene.reaction.association<-gsub("([[:alnum:]]+)\\.[[:digit:]]+"," \\1 ",RECON$Gene.reaction.association)
RECON$Gene.reaction.association <- sapply(RECON$Gene.reaction.association, function(gpr){
  woSpaces <- gsub("\\(|\\)|[[:blank:]]+","",gpr)
  paste0(unique(unlist(lapply(lapply(unlist(strsplit(woSpaces,"or")), function(gpr){strsplit(gpr,"and")}),function(gpr){paste0("( ",paste0(sort(unlist(gpr)),collapse = " and ")," )")}))),collapse = " or ")
},USE.NAMES = FALSE)

# Extrayendo las reacciones asociadas a los genes en RECON
Reactions <- sapply(unique(Astrocyte_Genes$ENTREZ_GENE[!is.na(Astrocyte_Genes$EC)]),function(enzyme){RECON$Formula[grep(paste0("[[:blank:]]",enzyme,"[[:blank:]]"),RECON$Gene.reaction.association)]})
Reactions <- unique(unlist(Reactions))

# GapFind y GapFill
Reactions <- gapFill(Reactions,RECON$Formula[RECON$Gene.reaction.association==""],consensus = TRUE)

# Añadiendo flujo
colnames(RECON) <- c("ID","DESCRIPTION","REACTION","GPR","REVERSIBLE","LOWER.BOUND","UPPER.BOUND","OBJECTIVE")
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

DMEM <- addReact(DMEM, id="MC", met=c("atp[c]","atp[e]"),
                 Scoef=c(-1,2), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("cys_L[e]","gthrd[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Reactions_Flux <- unique(c(Reactions_Flux,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

# Construyendo la reconstrucción
Astrocyte_Draft<- mapReactions(reactionList = unique(c(Reactions,Reactions_Flux)),referenceData = RECON,by = "REACTION")
convert2sbml(Astrocyte_Draft,"Results/Astrocyte_Draft.xml")

#
Astrocyte_DraftM <- readSBMLmod("Results/Astrocyte_Draft.xml")
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM)%in%react_id(findExchReact(Astrocyte_DraftM))] <- 0
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_ca2(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_glc(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_fe3(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_k(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_na1(e)'] <- -1000
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_HC02172(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_Rtotal(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_Rtotal2(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_Rtotal3(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_pi(e)'] <- -100
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_ala_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_arg_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_asn_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_asp_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_cys_L(e)'] <- -1000
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_glu_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_gly(e)'] <- -1000
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_his_L(e)'] <- -1000
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_ile_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_leu_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_lys_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_met_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_phe_L(e)'] <- -1000
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_pro_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_ser_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_thr_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_trp_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_val_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_btn(e)'] <- -100
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_chol(e)'] <- -1000
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_aqcobal(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_fol(e)'] <- -1000
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_inost(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_lnlc(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_lnlnca(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_ncam(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_pydxn(e)'] <- -1000
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_ribflv(e)'] <- -1000
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_thm(e)'] <- -1000
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_pyr(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_thymd(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_tyr_L(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_o2(e)'] <- -0.530
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_h2o(e)'] <- -100
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_cl(e)'] <- -1000
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_co2(e)'] <- 0.515
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_so4(e)'] <- -100
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_hdca(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_estradiol(e)'] <- -1
lowbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_nh4(e)'] <- -100
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM)%in%react_id(findExchReact(Astrocyte_DraftM))] <- 0
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_prostgd2(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_prostge1(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_prostge2(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_prostgf2(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_prostgh2(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_prostgi2(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_HC02202(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_o2(e)'] <- -0.515
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_HC02203(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_HC02204(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_HC02205(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_HC02206(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_HC02207(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_HC02208(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_HC02210(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_HC02213(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_HC02214(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_HC02216(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_HC02217(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_leuktrA4(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_leuktrB4(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_leuktrC4(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_leuktrD4(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_leuktrE4(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_leuktrF4(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_gthrd(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_co2(e)'] <- 0.530
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_lac_L(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_glc_D(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_gln_L(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_ser_D(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_atp(e)'] <- 1000
uppbnd(Astrocyte_DraftM)[react_id(Astrocyte_DraftM) == 'EX_taur(e)'] <- 1000

# 
lockedReactions <- function(model){
  if(!is.loaded("sybil")){require("sybil")}
  locked <- NULL
  pb <- txtProgressBar(min = 1,max = model@react_num,style=3)
  for (reaction in 1:model@react_num) {
    setTxtProgressBar(pb, reaction)
    model@obj_coef <- rep(0, model@react_num)
    model@obj_coef[reaction] <- 1
    FBA <- optimizeProb(model)
    locked <- unique(c(locked, model@react_id[as.vector(FBA@fluxdist@fluxes!=0)]))
  }
  close(pb)
  locked <- model@react_id[!model@react_id%in%locked]
  return(locked)
}
woFlux <- lockedReactions(Astrocyte_DraftM)

# 
Astrocyte_Draft$LOWER.BOUND <- Astrocyte_DraftM@lowbnd
Astrocyte_Draft$UPPER.BOUND <- Astrocyte_DraftM@uppbnd
Astrocyte_Reconstruction <- Astrocyte_Draft[!Astrocyte_Draft$ID%in%woFlux,]

#
convert2sbml(Astrocyte_Reconstruction,"Results/Astrocyte.xml")
