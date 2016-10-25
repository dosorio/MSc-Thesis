library(sybilSBML)

# Reading Model
matureAstrocyte_Model <- readSBMLmod("Results/matureAstrocyte.xml")

# Minimizing the total absolute fluxes MTF (Evaluando las múltiples posibles soluciones)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model)))

# Minimizing the total absolute fluxes MTF (Evaluando las múltiples posibles soluciones)
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_hdca(e)'] <- -0.209
uppbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_hdca(e)'] <- -0.209
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model)))


# Extraigo rutas metabolicas involucradas en la diferenciación de astrocito maduro
# AstrocyteMetabolism<-enrichPathway(unique(na.omit(Astrocyte_Genes$ENTREZ_GENE[Astrocyte_Genes$EC!=""])),organism = "human",readable = TRUE,pvalueCutoff = 0.05)
# ResumenAstrocyte<-summary(AstrocyteMetabolism)
