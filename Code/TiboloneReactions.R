# Tibolone Mecanism
lowbnd(DMEM)[react_id(DMEM) == 'EX_estradiol(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_crtsl(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_aldstrn(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_prgstrn(e)'] <- -1
lowbnd(DMEM)[react_id(DMEM) == 'EX_tststerone(e)'] <- -1

DMEM <- addReact(DMEM, id="MC", met=c("estradiol[e]","h2o2[m]","h2o[m]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- RECON[getFluxDist(optimizeProb(DMEM))!=0,3]

DMEM <- addReact(DMEM, id="MC", met=c("estradiol[e]","o2s[m]","h2o[m]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("crtsl[e]","h2o2[m]","h2o[m]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("crtsl[e]","o2s[m]","h2o[m]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("aldstrn[e]","o2s[m]","h2o[m]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("aldstrn[e]","h2o2[m]","h2o[m]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("prgstrn[e]","o2s[m]","h2o[m]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("prgstrn[e]","h2o2[m]","h2o[m]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("tststerone[e]","o2s[m]","h2o[m]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("tststerone[e]","h2o2[m]","h2o[m]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
DMEM <- addReact(DMEM, id="MC", met=c("estradiol[e]","h2o2[c]","h2o[c]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- RECON[getFluxDist(optimizeProb(DMEM))!=0,3]

DMEM <- addReact(DMEM, id="MC", met=c("estradiol[e]","o2s[c]","h2o[c]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("crtsl[e]","h2o2[c]","h2o[c]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("crtsl[e]","o2s[c]","h2o[c]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("aldstrn[e]","o2s[c]","h2o[c]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("aldstrn[e]","h2o2[c]","h2o[c]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("prgstrn[e]","o2s[c]","h2o[c]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("prgstrn[e]","h2o2[c]","h2o[c]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("tststerone[e]","o2s[c]","h2o[c]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))

DMEM <- addReact(DMEM, id="MC", met=c("tststerone[e]","h2o2[c]","h2o[c]"),
                 Scoef=c(-1,-1,1), reversible=FALSE,
                 lb=0, ub=1000, obj=1)
Tibolone <- unique(c(Tibolone,(RECON[getFluxDist(optimizeProb(DMEM))!=0,3])))
Tibolone <- mapReactions(reactionList = Tibolone[!Tibolone%in%Astrocyte_Reconstruction$REACTION],referenceData = RECON,by = "REACTION")
write.csv2(x = Tibolone,file = "Results/TiboloneReactions.csv",row.names = FALSE)