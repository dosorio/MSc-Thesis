require(sybilSBML)
require(Biobase)
require(exp2flux)
require(UniProt.ws)

# 
Astrocyte_Expression <- read.csv(file = "~/Documents/masterThesis/Data/GSE73721/GSE73721_Human_and_mouse_table.csv.gz",
                                 header = TRUE,
                                 row.names = "Gene")[,c(15:26)]
Astrocyte_Expression <- matrix(c(toupper(rownames(Astrocyte_Expression)),round(rowMeans(Astrocyte_Expression),3)),ncol = 2,dimnames = list(c(),c("GENECARDS","FPKM")))

# Extraigo ID's genes en diferentes bases de datos
Human <- (UniProt.ws(taxId=9606))
Astrocyte_Genes <- select(Human,keys=toupper(Astrocyte_Expression[,1]),columns = c("ENTREZ_GENE"),keytype ="GENECARDS")

# Merge
Astrocyte_ExpressionSet <- merge(Astrocyte_Genes,Astrocyte_Expression,by = "GENECARDS")
Astrocyte_ExpressionSet <- Astrocyte_ExpressionSet[complete.cases(Astrocyte_ExpressionSet),]
Astrocyte_ExpressionSet$FPKM <- sapply(Astrocyte_ExpressionSet$ENTREZ_GENE, function(ID){
  positions <- grep(paste0("^",ID,"$"),Astrocyte_ExpressionSet$ENTREZ_GENE)
  return(max(as.vector((Astrocyte_ExpressionSet$FPKM[positions]))))
})
Astrocyte_ExpressionSet <- unique(Astrocyte_ExpressionSet[,2:3])
Astrocyte_ExpressionSet <- matrix(as.numeric(Astrocyte_ExpressionSet$FPKM),dimnames = list(Astrocyte_ExpressionSet$ENTREZ_GENE,c("MatureAstrocytes")))
Astrocyte_ExpressionSet <- ExpressionSet(Astrocyte_ExpressionSet)


Astrocyte_Model <- readSBMLmod("Results/Astrocyte.xml") 
Astrocyte_Model

Astrocyte_Model <- exp2flux(Astrocyte_Model,Astrocyte_ExpressionSet,scale = TRUE)

plot(robAna(Astrocyte_Model,"EX_hdca(e)"))

# FBA
model_FBA <- optimizeProb(Astrocyte_Model)
model_FBA

# Minimizing the total absolute fluxes MTF (Evaluando las múltiples posibles soluciones)
model_MTF <- optimizeProb(Astrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
getNetFlux(getFluxDist(model_MTF,findExchReact(Astrocyte_Model)))

# Extraigo rutas metabolicas involucradas en la diferenciación de astrocito maduro
# AstrocyteMetabolism<-enrichPathway(unique(na.omit(Astrocyte_Genes$ENTREZ_GENE[Astrocyte_Genes$EC!=""])),organism = "human",readable = TRUE,pvalueCutoff = 0.05)
# ResumenAstrocyte<-summary(AstrocyteMetabolism)
