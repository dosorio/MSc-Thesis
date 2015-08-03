PCA_R<-princomp(~ TIME + CONFIDENT+ CONFPROT ,data = PCA,cor=TRUE)
biplot(PCA_R)
sort(PCA_R$scores[,1])
