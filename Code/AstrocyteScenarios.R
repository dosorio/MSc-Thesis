require(sybilSBML)
require(Biobase)
require(exp2flux)
require(UniProt.ws)

# Loading Gene Expression data
Astrocyte_Expression <- as.matrix(read.csv(file = "~/Documents/masterThesis/Data/GSE73721/GSE73721_Human_and_mouse_table.csv.gz",
                                 header = TRUE,row.names = "Gene")[,c(9:26)])

# Extraigo ID's genes en diferentes bases de datos
Human <- (UniProt.ws(taxId=9606))
Astrocyte_Genes <- select(Human,keys=rownames(Astrocyte_Expression),columns = c("ENTREZ_GENE"),keytype ="GENECARDS")

# Data Summary
for(ID in unique(Astrocyte_Genes$ENTREZ_GENE)){
  P <- Astrocyte_Genes$GENECARDS[Astrocyte_Genes$ENTREZ_GENE%in%ID]
  P <- P[P%in%rownames(Astrocyte_Expression)]
  if(length(P)>1){
    for (C in seq_len(ncol(Astrocyte_Expression))){
      Astrocyte_Expression[P,C] <- max(Astrocyte_Expression[P,C])
    }
    rownames(Astrocyte_Expression)[rownames(Astrocyte_Expression)%in%P] <- ID 
  }
  if(length(P)==1){
    rownames(Astrocyte_Expression)[rownames(Astrocyte_Expression)%in%P] <- ID 
  }
}
Astrocyte_Expression <- unique(Astrocyte_Expression)
Astrocyte_Expression <- Astrocyte_Expression[complete.cases(Astrocyte_Expression),]
Astrocyte_Expression <- Astrocyte_Expression[!is.na(rownames(Astrocyte_Expression)),]

# Creating healthy states
Astrocyte_Model <- readSBMLmod("Results/Astrocyte.xml")
Astrocyte_Model

par(mfcol=c(1,4))
# Fetal
Astrocyte_ExpressionSet <- ExpressionSet(Astrocyte_Expression[,c(1:6)])
fetalAstrocyte_Model <- exp2flux(Astrocyte_Model,Astrocyte_ExpressionSet,scale = FALSE)
plot(fluxVar(fetalAstrocyte_Model),ylim=c(-250,250),main=paste0("Fetal: ",round(optimizeProb(fetalAstrocyte_Model)@lp_obj,3)),ylab="Biomass Flux",xlab="Reaction")

# Young
Astrocyte_ExpressionSet <- ExpressionSet(Astrocyte_Expression[,c(7:9)])
youngAstrocyte_Model <- exp2flux(Astrocyte_Model,Astrocyte_ExpressionSet,scale = FALSE)
plot(fluxVar(youngAstrocyte_Model),ylim=c(-250,250),main=paste0("Young (8-16 yo): ",round(optimizeProb(youngAstrocyte_Model)@lp_obj,3)),ylab="Biomass Flux",xlab="Reaction")

# Adult
Astrocyte_ExpressionSet <- ExpressionSet(Astrocyte_Expression[,c(10:12)])
adultAstrocyte_Model <- exp2flux(Astrocyte_Model,Astrocyte_ExpressionSet,scale = FALSE)
plot(fluxVar(adultAstrocyte_Model),ylim=c(-250,250),main=paste0("Adult (21-35 yo): ",round(optimizeProb(adultAstrocyte_Model)@lp_obj,3)),ylab="Biomass Flux",xlab="Reaction")

# Mature
Astrocyte_ExpressionSet <- ExpressionSet(Astrocyte_Expression[,c(13:18)])
matureAstrocyte_Model <- exp2flux(Astrocyte_Model,Astrocyte_ExpressionSet,scale = FALSE)
plot(fluxVar(matureAstrocyte_Model),ylim=c(-250,250),main=paste0("Mature (47-63 yo): ",round(optimizeProb(matureAstrocyte_Model)@lp_obj,3)),ylab="Biomass Flux",xlab="Reaction")

# FBA
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA

# Minimizing the total absolute fluxes MTF (Evaluando las múltiples posibles soluciones)
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model)))

# Extraigo rutas metabolicas involucradas en la diferenciación de astrocito maduro
# AstrocyteMetabolism<-enrichPathway(unique(na.omit(Astrocyte_Genes$ENTREZ_GENE[Astrocyte_Genes$EC!=""])),organism = "human",readable = TRUE,pvalueCutoff = 0.05)
# ResumenAstrocyte<-summary(AstrocyteMetabolism)
