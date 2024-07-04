#Meth.TSS1500.aggr <- aggregate(x = Meth.TSS1500[,2:195], by = list(Meth.TSS1500$GeneSymbol), FUN = median)
library("samr")
rownames(Meth.TSS1500.aggr) <- Meth.TSS1500.aggr$Group.1
Meth.TSS1500.aggr$Group.1 <- NULL
colnames(Meth.TSS1500.aggr) <- gsub("\\.", "-", colnames(Meth.TSS1500.aggr))
sampleIDs <- colnames(Meth.TSS1500.aggr)
sampleIDs <- gsub("\\.", "-", sampleIDs)
samplesDat <- data.frame(matrix(nrow = length(sampleIDs), ncol = 7))
rownames(samplesDat) <- sampleIDs
for (j in 1:length(sampleIDs)) {
     tmpRow <- unlist(strsplit(sampleIDs[j], split = "-"))
     samplesDat[sampleIDs[j], ] <- tmpRow
}
sampleIDs1 <- as.character(samplesDat[, 4])
sampleIDs1 <- substr(sampleIDs1, 1, nchar(sampleIDs1) - 1)
sampleIDs1 <- as.numeric(sampleIDs1)
normalSamples <- rownames(samplesDat)[sampleIDs1 < 14 & sampleIDs1 > 9]
tumorSamples <- rownames(samplesDat)[sampleIDs1 < 3]
MethMat <- Meth.TSS1500.aggr[, c(normalSamples, tumorSamples)]
y <- c(rep(1,length(normalSamples)),rep(2,length(tumorSamples)))
data <- list(x=MethMat,y=y,logged2=FALSE,genenames=paste(row.names(MethMat)), geneid=paste(row.names(MethMat)))
samr.obj <- samr(data, resp.type = "Two class unpaired", nperms = 10000, random.seed = 123456)
delta.table <- samr.compute.delta.table(samr.obj, nvals=50)
del=1.5
siggenes.table <- samr.compute.siggenes.table(samr.obj, del, data, delta.table, min.foldchange = 2.0, all.genes=FALSE)
up <- siggenes.table$genes.up
down <- siggenes.table$genes.lo
upDown <- rbind(up,down)
low <- upDown[as.numeric(upDown[,8])<1,]
lowOrder <- low[order(low[,7], decreasing=TRUE),]
write.table(lowOrder, file = "MethTSS1500-Res.csv", sep = ",")
