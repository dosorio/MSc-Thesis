library(sybilSBML)

###########################################################################################
# HEALTHY
###########################################################################################

# Reading Model
matureAstrocyte_Model <- readSBMLmod("Results/matureAstrocyte.xml")

# 
dir.create("Results/healthyAstrocyte")
# BO Healthy
sink("Results/healthyAstrocyte/BO.txt")
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Glu-Gln
sink("Results/healthyAstrocyte/Glu-Gln.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("glu_L[e]","gln_L[e]"),
                 Scoef=c(-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model)))))
sink()

# Gly-serD
sink("Results/healthyAstrocyte/Gly-SerD.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("gly[e]","ser_D[e]"),
                                  Scoef=c(-1,1), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Glc-Lactate
sink("Results/healthyAstrocyte/Glc-Lactate.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("glc_D[e]","lac_L[e]"),
                                  Scoef=c(-1,2), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()


# Glc-ATP
sink("Results/healthyAstrocyte/Glc-ATP.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("glc_D[e]","atp[e]"),
                                  Scoef=c(-1,36), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Cys-GTHRD
sink("Results/healthyAstrocyte/Cys-GTHRD.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("cys_L[e]","glu_L[c]","gly[c]","gthrd[e]"),
                                  Scoef=c(-1,-1,-1,1), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

###########################################################################################
# INFLAMMATED 
###########################################################################################

# Reading Model
matureAstrocyte_Model <- readSBMLmod("Results/matureAstrocyte.xml")

# Adding Palmitate
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_hdca(e)'] <- -0.226
uppbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_hdca(e)'] <- -0.226


# Minimizing the total absolute fluxes MTF (Evaluando las múltiples posibles soluciones)
dir.create("Results/inflammatedAstrocyte")

# BO Inflammated
sink("Results/inflammatedAstrocyte/BO.txt")
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Glu-Gln
sink("Results/inflammatedAstrocyte/Glu-Gln.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("glu_L[e]","gln_L[e]"),
                                  Scoef=c(-1,1), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Gly-serD
sink("Results/inflammatedAstrocyte/Gly-SerD.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("gly[e]","ser_D[e]"),
                                  Scoef=c(-1,1), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Glc-Lactate
sink("Results/inflammatedAstrocyte/Glc-Lactate.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("glc_D[e]","lac_L[e]"),
                                  Scoef=c(-1,2), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()


# Glc-ATP
sink("Results/inflammatedAstrocyte/Glc-ATP.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("glc_D[e]","atp[e]"),
                                  Scoef=c(-1,36), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Cys-GTHRD
sink("Results/inflammatedAstrocyte/Cys-GTHRD.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("cys_L[e]","glu_L[c]","gly[c]","gthrd[e]"),
                                  Scoef=c(-1,-1,-1,1), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

###########################################################################################
# INFLAMMATED + TIBOLONE
###########################################################################################

# Reading Model
matureAstrocyte_Model <- readSBMLmod("Results/matureAstrocyte.xml")

# Adding Palmitate
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_hdca(e)'] <- -0.226
uppbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_hdca(e)'] <- -0.226
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_estradiol(e)'] <- -1
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_prgstrn(e)'] <- -1
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_tststerone(e)'] <- -1
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_crtsl(e)'] <- -1
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_aldstrn(e)'] <- -1

# Adding Tibolone
matureAstrocyte_Model <- addExchReact(matureAstrocyte_Model,met = "tib[e]",lb = -1,ub = -1)
matureAstrocyte_Model <- addReact(model = matureAstrocyte_Model,
                                  id = "TiboloneMetabolism",
                                  met = c("tib[e]","OHtib_3a[e]","OHtib_3b[e]","d4_tib[e]"),
                                  Scoef = c(-1,1,1,1),
                                  reversible = FALSE,
                                  ub = 1000)
matureAstrocyte_Model <- addReact(model = matureAstrocyte_Model,
                                  id = "EstrogenAgonist1",
                                  met = c("OHtib_3a[e]","estradiol[c]"),
                                  Scoef = c(-1,1),
                                  reversible = FALSE,
                                  ub = 1000)
matureAstrocyte_Model <- addReact(model = matureAstrocyte_Model,
                                  id = "EstrogenAgonist2",
                                  met = c("OHtib_3b[e]","estradiol[c]"),
                                  Scoef = c(-1,1),
                                  reversible = FALSE,
                                  ub = 1000)
matureAstrocyte_Model <- addReact(model = matureAstrocyte_Model,
                                  id = "ProgesteroneAndAndrogenicAgonist",
                                  met = c("d4_tib[e]","prgstrn[e]","tststerone[e]","prgstrn[c]","tststerone[c]"),
                                  Scoef = c(-1,-2,-1,2,1),
                                  reversible = FALSE,
                                  ub = 1000)
matureAstrocyte_Model <- addReact(model = matureAstrocyte_Model,
                                  id = "ProgesteroneAndAndrogenicAntagonist",
                                  met = c("OHtib_3a[e]","OHtib_3b[e]","prgstrn[c]","tststerone[c]","prgstrn[e]","tststerone[e]"),
                                  Scoef = c(-1,-1,-1,-1,1,1),
                                  reversible = FALSE,
                                  ub = 1000)
matureAstrocyte_Model <- addReact(model = matureAstrocyte_Model,
                                  id = "MineralocorticoidAntagonist",
                                  met = c("OHtib_3a[e]","OHtib_3b[e]","crtsl[c]","aldstrn[c]","crtsl[e]","aldstrn[e]"),
                                  Scoef = c(-1,-1,-1,-1,1,1),
                                  reversible = FALSE,
                                  ub = 1000)

plot(robAna(matureAstrocyte_Model,"Ex_tib[e]",rng = c(-1,0)))
# Minimizing the total absolute fluxes MTF (Evaluando las múltiples posibles soluciones)
dir.create("Results/inflammatedAstrocyte")

# BO Inflammated
sink("Results/inflammatedAstrocyte/BO.txt")
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Glu-Gln
sink("Results/inflammatedAstrocyte/Glu-Gln.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("glu_L[e]","gln_L[e]"),
                                  Scoef=c(-1,1), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Gly-serD
sink("Results/inflammatedAstrocyte/Gly-SerD.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("gly[e]","ser_D[e]"),
                                  Scoef=c(-1,1), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Glc-Lactate
sink("Results/inflammatedAstrocyte/Glc-Lactate.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("glc_D[e]","lac_L[e]"),
                                  Scoef=c(-1,2), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()


# Glc-ATP
sink("Results/inflammatedAstrocyte/Glc-ATP.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("glc_D[e]","atp[e]"),
                                  Scoef=c(-1,36), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Cys-GTHRD
sink("Results/inflammatedAstrocyte/Cys-GTHRD.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("cys_L[e]","glu_L[c]","gly[c]","gthrd[e]"),
                                  Scoef=c(-1,-1,-1,1), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()


# Extraigo rutas metabolicas involucradas en la diferenciación de astrocito maduro
# AstrocyteMetabolism<-enrichPathway(unique(na.omit(Astrocyte_Genes$ENTREZ_GENE[Astrocyte_Genes$EC!=""])),organism = "human",readable = TRUE,pvalueCutoff = 0.05)
# ResumenAstrocyte<-summary(AstrocyteMetabolism)
