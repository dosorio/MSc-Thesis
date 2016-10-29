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
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_hdca(e)'] <- -0.208
uppbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_hdca(e)'] <- -0.208


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
matureAstrocyte_Model <- read.csv2("Results/matureAstrocyte.csv")
Tibolone <- read.csv2("Results/matureTiboloneReactions.csv")
matureAstrocyte_Model <- rbind(matureAstrocyte_Model,Tibolone)
convert2sbml(matureAstrocyte_Model,"Results/matureAstrocyte_Tibolone.xml")
matureAstrocyte_Model <- readSBMLmod("Results/matureAstrocyte_Tibolone.xml")

# Adding Palmitate
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_hdca(e)'] <- -0.208
uppbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_hdca(e)'] <- -0.208
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_estradiol(e)'] <- -1
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_prgstrn(e)'] <- -1
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_tststerone(e)'] <- -1
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_crtsl(e)'] <- -1
lowbnd(matureAstrocyte_Model)[react_id(matureAstrocyte_Model) == 'EX_aldstrn(e)'] <- -1

dir.create("Results/Tibolone")

# BO Inflammated
sink("Results/Tibolone/BO.txt")
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Essential reactions
oV <- optimizeProb(matureAstrocyte_Model)@lp_obj
E <- NULL
for(reactionID in Tibolone$ID){
  lb <- lowbnd(matureAstrocyte_Model)[matureAstrocyte_Model@react_id==reactionID]
  ub <- uppbnd(matureAstrocyte_Model)[matureAstrocyte_Model@react_id==reactionID]
  lowbnd(matureAstrocyte_Model)[matureAstrocyte_Model@react_id==reactionID] <- 0
  uppbnd(matureAstrocyte_Model)[matureAstrocyte_Model@react_id==reactionID] <- 0
  if(optimizeProb(matureAstrocyte_Model)@lp_obj < (oV - (oV*0.1))){
    E <- c(E,reactionID)
  }
  lowbnd(matureAstrocyte_Model)[matureAstrocyte_Model@react_id==reactionID] <- 0
  uppbnd(matureAstrocyte_Model)[matureAstrocyte_Model@react_id==reactionID] <- 0
}

# Glu-Gln
sink("Results/Tibolone/Glu-Gln.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("glu_L[e]","gln_L[e]"),
                                  Scoef=c(-1,1), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Gly-serD
sink("Results/Tibolone/Gly-SerD.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("gly[e]","ser_D[e]"),
                                  Scoef=c(-1,1), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Glc-Lactate
sink("Results/Tibolone/Glc-Lactate.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("glc_D[e]","lac_L[e]"),
                                  Scoef=c(-1,2), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()


# Glc-ATP
sink("Results/Tibolone/Glc-ATP.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("glc_D[e]","atp[e]"),
                                  Scoef=c(-1,36), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

# Cys-GTHRD
sink("Results/Tibolone/Cys-GTHRD.txt")
matureAstrocyte_Model <- addReact(matureAstrocyte_Model, id="MC", met=c("cys_L[e]","glu_L[c]","gly[c]","gthrd[e]"),
                                  Scoef=c(-1,-1,-1,1), reversible=FALSE,
                                  lb=0, ub=1000, obj=1)
model_FBA <- optimizeProb(matureAstrocyte_Model)
model_FBA
model_MTF <- optimizeProb(matureAstrocyte_Model, algorithm = "mtf", wtobj = mod_obj(model_FBA))
print(getNetFlux(getFluxDist(model_MTF,findExchReact(matureAstrocyte_Model))))
sink()

#### Tibolone

# Extraigo rutas metabolicas involucradas en la diferenciación de astrocito maduro
# AstrocyteMetabolism<-enrichPathway(unique(na.omit(Astrocyte_Genes$ENTREZ_GENE[Astrocyte_Genes$EC!=""])),organism = "human",readable = TRUE,pvalueCutoff = 0.05)
# ResumenAstrocyte<-summary(AstrocyteMetabolism)
