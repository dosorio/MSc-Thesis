library(exp2flux)
library(gage)
library(gridExtra)
library(minval)
library(sybilSBML)
library(lattice)

Astrocyte <- read.csv2("Results/Astrocyte.csv")
length(grep("EX_",Astrocyte$ID))
table(Astrocyte$GPR[lengths(sapply(Astrocyte$REACTION[!grepl("EX_",Astrocyte$ID)], compartments))>1]!="")
RXN <- AGenes[AGenes$ENTREZ_GENE%in%unique(unlist(strsplit(gsub("\\(|or|and|\\)","",Astrocyte$GPR[lengths(sapply(Astrocyte$REACTION[!grepl("EX_",Astrocyte$ID)], compartments))==1]),"[[:blank:]]+"))),]

#pdf("Documents/thesisDocument/neuroprotective/RXN.pdf",width = 8,height = 8)
par(mfcol=c(1,3))
pie(c(60,1080,1607),
    labels = c("2.2%\nExchange\nReactions",
               "39.3%\nTransport Reactions\n",
               "\n\n58.5%\nCompartmentalized\nReactions"),main="A: Type of Reactions")

pie(c(407,138,291,126,54,17,28),
    labels=c("25.3%\nExpontaneous\nReactions\n\n",
             "15.8%\nOxido-\nreductases",
             "33.2%\nTransferases",
             "14.4%\nHydrolases",
             "6.2%\nLyases",
             "1.9% Isomerases",
             "2.3% Ligases"),main="B: Catalytic Activity Required")

pie(table(sapply(Astrocyte$REACTION[lengths(sapply(Astrocyte$REACTION[!grepl("EX_",Astrocyte$REACTION)],compartments))==1],compartments)),
    labels = c("39.1%\nCytoplasm\n",
               "5.5%\nExtracellular",
               "8.7%\nGolgi\napparatus",
               "6.3%\nLysosome",
               "\n20.9%\nMitochondrion",
               "\n4.2%\nNucleus",
               "\n6.9%\nEndoplasmic\nReticulum",
               "8.5%\nPeroxisome"),main="C: Reactions Distribution by Compartment")

dev.off()
pdf("Documents/thesisDocument/neuroprotective/Pathways.pdf",height = 15,width = 15)
par(las=2,mar=c(4,30,5,5),cex=0.9)
# data <- try(kegg.gsets(species = "hsa", id.type = "entrez"))
# data <- matrix(gsub("[[:digit:]]+$","",names(unlist(data$kg.sets))),dimnames = list(as.vector(unlist(data$kg.sets)),c()))
# data[,1] <- gsub("hsa[[:digit:]]+ ","",data[,1])
# genes <- unique(unlist(strsplit(gsub("\\(|or|and|\\)","",Astrocyte$GPR)," ")))
d <- sort(table(data[as.character(genes[genes%in%rownames(data)]),])/sum(table(data[as.character(genes[genes%in%rownames(data)]),]))*100,decreasing = FALSE)
barplot(d,horiz = TRUE,main = "Percentage of Reactions by Metabolic Pathway",xlab="% Reactions")
dev.off()

# par(mfcol=c(1,2))
healthy <- readSBMLmod("Results/matureAstrocyte.xml")
model_FBA <- optimizeProb(healty)
model_MTF <- optimizeProb(healty, algorithm = "mtf", wtobj = mod_obj(model_FBA))
h <- getNetFlux(getFluxDist(model_MTF,findExchReact(healty)))
inflammated <- healthy
lowbnd(inflammated)[inflammated@react_id=="EX_hdca(e)"] <- -0.208
uppbnd(inflammated)[inflammated@react_id=="EX_hdca(e)"] <- -0.208
model_FBA <- optimizeProb(inflammated)
model_MTF <- optimizeProb(inflammated, algorithm = "mtf", wtobj = mod_obj(model_FBA))
i <- getNetFlux(getFluxDist(model_MTF,findExchReact(inflammated)))
medicated <- readSBMLmod("Results/matureAstrocyte_Tibolone.xml")
lowbnd(medicated)[medicated@react_id=="EX_hdca(e)"] <- -0.208
uppbnd(medicated)[medicated@react_id=="EX_hdca(e)"] <- -0.208
model_FBA <- optimizeProb(medicated)
model_MTF <- optimizeProb(medicated, algorithm = "mtf", wtobj = mod_obj(model_FBA))
m <- getNetFlux(getFluxDist(model_MTF,findExchReact(medicated)))
h <- cbind(h@react_id,h@rate)
i <- cbind(i@react_id,i@rate)
m <- cbind(m@react_id,m@rate)
inputs <- merge(merge(m,i,by="V1"),h,by="V1")
inputs <- as.matrix(inputs[,2:4])
rownames(inputs) <- c("L-alanine","L-arginine","L-asparagine","L-aspartate","ATP","Calcium","Choline","Chloride",
                      "Carbon Dioxide","L-cysteine","Fe2+","Fe3+","Folate","D-glucose","L-glutamine","L-glutamate",
                      "Glycine","Reduced Glutathione","Water","Prostaglandin-D1","Prostaglandin-D3","Prostaglandin-E3","Prostaglandin-F1a",
                      "Prostaglandin-F2b","Prostaglandin-G2","Palmitate","L-histidine","L-isoleucine","Myo-inositol","(S)-lactate",
                      "L-leucine","Leukotriene-B4","Leukotriene-E4","Linoleate","aLinolenate","L-lysine","L-methionine",
                      "Sodium","Nicotinamide","Ammonium","Nitric Oxide","Oxygen","L-phenylalanine","Hydrogenphosphate","L-proline",
                      "Prostaglandin-D2","Prostaglandin-E1","Prostaglandin-H2","Prostaglandin-I2","Riboflavin","RECON2-Rtotal",
                      "RECON2-Rtotal2","D-serine","L-serine","Sulfate","L-threonine","Thymidine","L-tryptophan",
                      "L-tyrosine","L-valine")

colnames(inputs) <- c("Medicated","Inflammated","Healthy")
levelplot((inputs),at=seq(-3,3, length.out=96),colorkey=list(space="top",labels=list(at=c(-3,-2,-1,0,1,2,3),labels=c("Uptake\n-3","-2","-1","0","1","2","Release\n3"))), scales=list(x=list(rot=90)),ylab="Scenario",xlab="Exchange Metabolites")
# plot(ifV,col=ifelse(C==TRUE,"black","red"),ylab="Flux",xlab="Reaction")
# medicated <- readSBMLmod("Results/matureAstrocyte_Tibolone.xml")
# lowbnd(medicated)[medicated@react_id=="EX_hdca(e)"] <- -0.208
# uppbnd(medicated)[medicated@react_id=="EX_hdca(e)"] <- -0.208
# mfV <- fluxVar(medicated)
# C <- (ifV@lp_obj!=mfV@lp_obj[1:healty@react_num])[1:healty@react_num] & (ifV@lp_obj!=mfV@lp_obj[(healty@react_num+1):(healty@react_num*2)])[(healty@react_num+1):(healty@react_num*2)]
# table(C)
# plot(ifV,col=ifelse(C==TRUE,"black","blue"),ylab="Flux",xlab="Reaction")


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
  genesD <- rownames(differences)[differences[,3]<0]
  genesD <- unique(unlist(model1@genes[model1@react_id%in%genesD]))
  genesD <- genesD[genesD%in%rownames(metabolicPathways)]
  genesD <- table(metabolicPathways[genesD,])
  genesK <- table(metabolicPathways)
  genesK <- genesK[names(genesK)%in%names(genesD)]
  genesD <- sort((genesD/genesK)*100)
  par(mfcol=c(1,2),las=1,mar=c(3,15,3,3))
  genesAD <- genesA[!names(genesA)%in%names(genesD)]
  genesDD <- genesD[!names(genesD)%in%names(genesA)]
  barplot(genesAD,horiz = TRUE,main = "Pathway % Activation",cex.names = 0.7,adj=1)
  barplot(genesDD,horiz = TRUE,main = "Pathway % Inactivation",cex.names = 0.7,adj=1)
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

