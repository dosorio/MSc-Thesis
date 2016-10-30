library(exp2flux)
library(gage)
library(gridExtra)

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
  barplot(genesA,horiz = TRUE,main = "Pathway % Activation",cex.names = 0.7)
  genesD <- rownames(differences)[differences[,3]<0]
  genesD <- unique(unlist(model1@genes[model1@react_id%in%genesD]))
  genesD <- genesD[genesD%in%rownames(metabolicPathways)]
  genesD <- table(metabolicPathways[genesD,])
  genesK <- table(metabolicPathways)
  genesK <- genesK[names(genesK)%in%names(genesD)]
  genesD <- sort((genesD/genesK)*100)
  barplot(genesD,horiz = TRUE,main = "Pathway % Inactivation",cex.names = 0.7)
}
pdf(file = "Slides/Figures/Healthy2Inflammated.pdf",width = 10,height = 6.5)
metabolicChanges(healthy,inflammated)
dev.off()
h2i <- fluxDifferences(healthy,inflammated)
i2t <- fluxDifferences(inflammated,tibolone)