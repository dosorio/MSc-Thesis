require(sybilSBML)
require(Biobase)
require(exp2flux)
require(UniProt.ws)
require(exp2flux)
require(minval)

# Loading Gene Expression data
Astrocyte_Expression <- as.matrix(read.csv(file = "Data/GSE73721/GSE73721_Human_and_mouse_table.csv.gz",
                                           header = TRUE,row.names = "Gene")[,c(9:26)])

# Extraigo ID's genes en diferentes bases de datos
Human <- (UniProt.ws(taxId=9606))
Astrocyte_Genes <- select(Human,keys=rownames(Astrocyte_Expression),columns = c("ENTREZ_GENE","EC"),keytype ="GENECARDS")

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
Astrocyte_Modelcsv <- read.csv2("Results/Astrocyte.csv")
Astrocyte_Model <- readSBMLmod("Results/Astrocyte.xml")
Astrocyte_Model

# División del gráfico
pdf(file = "Slides/Figures/Astrocyte_MetabolicChanges.pdf",width = 10,height = 6.5)
par(mfcol=c(1,4))

# Fetal
Astrocyte_ExpressionSet <- ExpressionSet(Astrocyte_Expression[,c(1:6)])
fetalAstrocyte_Model <- exp2flux(Astrocyte_Model,Astrocyte_ExpressionSet,scale = FALSE)
plot(fluxVar(fetalAstrocyte_Model),ylim=c(-250,250),main=paste0("Fetal: ",round(optimizeProb(fetalAstrocyte_Model)@lp_obj,3)),ylab="Biomass Flux",xlab="Reaction")
fetalAstrocyte_Modelcsv <- Astrocyte_Modelcsv
fetalAstrocyte_Modelcsv$LOWER.BOUND <- fetalAstrocyte_Model@lowbnd
fetalAstrocyte_Modelcsv$UPPER.BOUND <- fetalAstrocyte_Model@uppbnd
write.csv2(fetalAstrocyte_Modelcsv,file = "Results/fetalAstrocyte.csv",row.names = FALSE)
convert2sbml(fetalAstrocyte_Modelcsv,"Results/fetalAstrocyte.xml")

# Young
Astrocyte_ExpressionSet <- ExpressionSet(Astrocyte_Expression[,c(7:9)])
youngAstrocyte_Model <- exp2flux(Astrocyte_Model,Astrocyte_ExpressionSet,scale = FALSE)
plot(fluxVar(youngAstrocyte_Model),ylim=c(-250,250),main=paste0("Young (8-16 yo): ",round(optimizeProb(youngAstrocyte_Model)@lp_obj,3)),ylab="Biomass Flux",xlab="Reaction")
youngAstrocyte_Modelcsv <- Astrocyte_Modelcsv
youngAstrocyte_Modelcsv$LOWER.BOUND <- youngAstrocyte_Model@lowbnd
youngAstrocyte_Modelcsv$UPPER.BOUND <- youngAstrocyte_Model@uppbnd
write.csv2(youngAstrocyte_Modelcsv,file = "Results/youngAstrocyte.csv",row.names = FALSE)
convert2sbml(youngAstrocyte_Modelcsv,"Results/youngAstrocyte.xml")

# Adult
Astrocyte_ExpressionSet <- ExpressionSet(Astrocyte_Expression[,c(10:12)])
adultAstrocyte_Model <- exp2flux(Astrocyte_Model,Astrocyte_ExpressionSet,scale = FALSE)
plot(fluxVar(adultAstrocyte_Model),ylim=c(-250,250),main=paste0("Adult (21-35 yo): ",round(optimizeProb(adultAstrocyte_Model)@lp_obj,3)),ylab="Biomass Flux",xlab="Reaction")
adultAstrocyte_Modelcsv <- Astrocyte_Modelcsv
adultAstrocyte_Modelcsv$LOWER.BOUND <- adultAstrocyte_Model@lowbnd
adultAstrocyte_Modelcsv$UPPER.BOUND <- adultAstrocyte_Model@uppbnd
write.csv2(adultAstrocyte_Modelcsv,file = "Results/adultAstrocyte.csv",row.names = FALSE)
convert2sbml(adultAstrocyte_Modelcsv,"Results/adultAstrocyte.xml")

# Mature
Astrocyte_ExpressionSet <- ExpressionSet(Astrocyte_Expression[,c(13:18)])
matureAstrocyte_Model <- exp2flux(Astrocyte_Model,Astrocyte_ExpressionSet,scale = FALSE)
plot(fluxVar(matureAstrocyte_Model),ylim=c(-250,250),main=paste0("Mature (47-63 yo): ",round(optimizeProb(matureAstrocyte_Model)@lp_obj,3)),ylab="Biomass Flux",xlab="Reaction")
matureAstrocyte_Modelcsv <- Astrocyte_Modelcsv
matureAstrocyte_Modelcsv$LOWER.BOUND <- matureAstrocyte_Model@lowbnd
matureAstrocyte_Modelcsv$UPPER.BOUND <- matureAstrocyte_Model@uppbnd
write.csv2(matureAstrocyte_Modelcsv,file = "Results/matureAstrocyte.csv",row.names = FALSE)
convert2sbml(matureAstrocyte_Modelcsv,"Results/matureAstrocyte.xml")

## Tibolone
Tibolone <- read.csv2("Results/TiboloneReactions.csv")
convert2sbml(Tibolone,"Results/Tibolone.xml")
Tibolone_Model <- readSBMLmod("Results/Tibolone.xml")
Tibolone_Model <- exp2flux(Tibolone_Model,Astrocyte_ExpressionSet)
Tibolone$LOWER.BOUND <- Tibolone_Model@lowbnd
Tibolone$UPPER.BOUND <- Tibolone_Model@uppbnd
write.csv2(Tibolone,file = "Results/matureTiboloneReactions.csv",row.names = FALSE)
dev.off()