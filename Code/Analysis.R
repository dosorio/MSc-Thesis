library(exp2flux)
library(gage)
library(gridExtra)
library(minval)
library(sybilSBML)

Astrocyte <- read.csv2("Results/Astrocyte.csv")
length(grep("EX_",Astrocyte$ID))
table(Astrocyte$GPR[lengths(sapply(Astrocyte$REACTION[!grepl("EX_",Astrocyte$ID)], compartments))>1]!="")
RXN <- AGenes[AGenes$ENTREZ_GENE%in%unique(unlist(strsplit(gsub("\\(|or|and|\\)","",Astrocyte$GPR[lengths(sapply(Astrocyte$REACTION[!grepl("EX_",Astrocyte$ID)], compartments))==1]),"[[:blank:]]+"))),]

metabolicChanges <- function(model1,model2,main){
  differences <- fluxDifferences(model1,model2)
  metabolicPathways <- kegg.gsets(species = "hsa", id.type = "kegg")
  metabolicPathways <- matrix(gsub("[[:digit:]]+$","",names(unlist(metabolicPathways$kg.sets))),dimnames = list(as.vector(unlist(metabolicPathways$kg.sets)),c()))
  metabolicPathways[,1] <- gsub("hsa[[:digit:]]+[[:blank:]]+","",metabolicPathways[,1])
  genesA <- rownames(differences)[differences[,3]>0]
  genesA <- unique(unlist(model1@genes[model1@react_id%in%genesA]))
  genesA <- genesA[genesA%in%rownames(metabolicPathways)]
  genesA <- table(metabolicPathways[genesA,])
  genesK <- table(metabolicPathways)
  genesK <- genesK[names(genesK)%in%names(genesA)]
  genesA <- sort((genesA/genesK)*100)
  par(mfcol=c(1,2),las=1,mar=c(3,15,3,3))
  barplot(genesA,horiz = TRUE,main = "Pathway % Activation",cex.names = 0.7,adj=1)
  genesD <- rownames(differences)[differences[,3]<0]
  genesD <- unique(unlist(model1@genes[model1@react_id%in%genesD]))
  genesD <- genesD[genesD%in%rownames(metabolicPathways)]
  genesD <- table(metabolicPathways[genesD,])
  genesK <- table(metabolicPathways)
  genesK <- genesK[names(genesK)%in%names(genesD)]
  genesD <- sort((genesD/genesK)*100)
  barplot(genesD,horiz = TRUE,main = "Pathway % Inactivation",cex.names = 0.7,adj=1)
}
pdf(file = "Slides/Figures/Healthy2Inflammated.pdf",width = 10,height = 6.5)
metabolicChanges(healthy,inflammated)
dev.off()
pdf(file = "Slides/Figures/Inflammated2Tibolone.pdf",width = 10,height = 6.5)
metabolicChanges(inflammated,tibolone)
dev.off()

h2i <- fluxDifferences(healthy,inflammated)
i2t <- fluxDifferences(inflammated,tibolone)
pdf(file = "Slides/Figures/Effects.pdf",width = 10,height = 6.5,encoding = 'ISOLatin2.enc')
par(mfcol=c(1,6),las=2,mar=c(7,3,4,3))
barplot(c(Healthy=0.377,Inflammated=0.318,Tibolone=0.427),col=c("darkblue","red","dodgerblue2"),main="BIOMASS\nDMEM medium")
barplot(c(Healthy=4.652,Inflammated=1.893,Tibolone=1.893),col=c("darkblue","red","dodgerblue2"),main="CYS[e] to GTHRD[e]\nConvertion")
barplot(c(Healthy=1.813,Inflammated=0.507,Tibolone=0.507),col=c("darkblue","red","dodgerblue2"),main="GLC[e] to ATP[e]\nConvertion")
barplot(c(Healthy=1.792,Inflammated=0.458,Tibolone=0.458),col=c("darkblue","red","dodgerblue2"),main="GLC[e] to LACTATE[e]\nConvertion")
barplot(c(Healthy=3.185,Inflammated=1.029,Tibolone=1.029),col=c("darkblue","red","dodgerblue2"),main="GLU[e] to GLN[e]\nConvertion")
barplot(c(Healthy=1000.377,Inflammated=1000.377,Tibolone=1000.377),col=c("darkblue","red","dodgerblue2"),main="GLY[e] to dSER[e]\nConvertion")
dev.off()

