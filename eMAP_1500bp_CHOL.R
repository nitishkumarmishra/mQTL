#library(edgeR)
library(edgeR)
library(eMap) # see references in paper for web page
library(quantsmooth)
load("BMIQ_Meth.Rda")
res <- read.csv("refseq_genes_hg19_finished.txt", header=TRUE, sep="\t")
probeinfo <- read.csv("Probeinfo.txt", header=TRUE, sep="\t")
geneExp <- read.csv("CHOL.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt", header = TRUE, sep = "\t", stringsAsFactors=FALSE, row.names = 1, check.names = FALSE)
geneExp <- geneExp[-c(1),which(geneExp[1,]=="raw_count")]
rownames(geneExp) <- gsub(pattern = "SLC35E2\\|728661", "SLC35E2B\\|728661", x = rownames(geneExp))
colnames(geneExp) <- substr(colnames(geneExp),1,16)
geneExp <- geneExp[-grep("\\?", rownames(geneExp)),]
geneExp <- round(data.matrix(geneExp),0)
cancerID <- grep("01A", colnames(geneExp))
normalID <- grep("11A", colnames(geneExp))
cancerExp <- geneExp[,cancerID]
cnts <- cancerExp[rowSums(cancerExp==0)< ncol(cancerExp)*0.2,] ## Remove all gene which have 25% zero's
keep <- rowSums(cpm(cnts)>1) >= ncol(cnts)*0.25 #### Select only genes which have have CPM > 1 for >=50% samples
cnts <- cnts[keep,]

cnts <- cbind(cnts[,cancerID], cnts[,normalID])
design <- model.matrix(~0 + factor(c(rep(2, length(cancerID)), rep(1, length(normalID)))))
colnames(design) <- c("Normal", "Tumor")
#cont.matrix <- makeContrasts("Tumor-Normal", levels = design)
#cnf <- calcNormFactors(cnts, method = "TMM")
v <- voom(cnts, design, plot = TRUE)
expr <- as.data.frame(v$E)
expr$gene <- sapply(strsplit(rownames(expr),"\\|"),'[[',1)

exprMerge <- merge(expr, res, by.x="gene", by.y="name2")
exprMerge$chrom <- gsub("chr","", exprMerge$chrom)
colnames(exprMerge)[48] <- c("Chr")
colnames(exprMerge)[50] <- c("TSS")
expr <- exprMerge[order(exprMerge$Chr,exprMerge$TSS),]
expr<-expr[order(expr$Chr,expr$TSS),]
me<-as.matrix(expr[,2:46]);rownames(me)<-expr[,1]
eChr <- expr$Chr
eChr <- gsub("X", "23", eChr, perl = TRUE)
eChr <- gsub("Y", "24", eChr, perl = TRUE)
eChr<-as.integer(eChr)
ePos<-as.numeric(expr$TSS)

meth<-cbind(IlmnID=rownames(BMIQ.Meth),data.frame(BMIQ.Meth)) # Where "beta" is a matrix of methylation values. 
meth<-merge(probeinfo,meth,by.x="Probe",by.y="IlmnID",all.x=F,all.y=F)
meth<-meth[order(meth$Chr,meth$Pos),]
mm<-as.matrix(meth[,-(1:7)]);rownames(mm)<-meth[,1]
colnames(mm) <- gsub("\\.", "-", colnames(mm))
mm<-mm[,colnames(mm)%in%colnames(me)]
me<-me[,colnames(me)%in%colnames(mm)]
me<-me[,order(colnames(me))]
mm<-mm[,order(colnames(mm))]
mChr<-meth$Chr
mChr <- gsub("X", "23", mChr, perl = TRUE)#mChr[mChr=="X"]<-23  #I used gsub
mChr <- gsub("Y", "24", mChr, perl = TRUE)#mChr[mChr=="Y"]<-24 
mChr<-as.integer(mChr)
mPos<-as.numeric(meth$Pos)
#save.image("CHOLeMAP.RData")

eMap1(me=me,mm1=mm,output.tag="output",p.cut=1,cis.only=T,cis.distance=1500,eChr=eChr,ePos=ePos,mChr=mChr,mPos=mPos)
eqtl<-read.table("output_eqtl.txt",header=T,sep="\t",na.strings="NA",row.names=NULL)
eqtl<-eqtl[,colnames(eqtl)%in%c("Gene_ID","Marker_ID","b1","b1_p")]
eqtl$bonf_p<-p.adjust(eqtl$b1_p,method="bonferroni")
eannot<-cbind(data.frame(nr=1:nrow(me)),Gene=rownames(me),eChr=eChr,ePos=ePos)
mannot<-cbind(data.frame(nr=1:nrow(mm)),Probe=rownames(mm),mChr=mChr,mPos=mPos)
eqtl<-merge(eqtl,eannot,all.x=T,all.y=F,by.x="Gene_ID",by.y="nr")
eqtl<-merge(eqtl,mannot,all.x=T,all.y=F,by.x="Marker_ID",by.y="nr")

eqtl <- merge(eqtl,expr[,colnames(expr)%in%c("gene","strand")],all.x=T,all.y=F,by.x="Gene",by.y="gene")
eqtl$dist<-(eqtl$mPos-eqtl$ePos)*eqtl$strand
saveRDS(eqtl,"output_1500bp_eqtl_annot.rds")
eqtl<-eqtl[eqtl$bonf_p<0.05,]
saveRDS(eqtl,"output_1500bp_eqtl_annot_bonf.rds")
write.csv(eqtl, file="eqtl_1500bp-0.5.csv")

############## Positive correlations ##############
pick<-eqtl$b1>0
logp<-(-log10(as.numeric(eqtl$b1_p[pick])))
mpos<-as.numeric(eqtl$mPos[pick])
tss<-as.numeric(eqtl$ePos[pick])
strand<-as.integer(eqtl$strand[pick])
disttss<-(mpos-tss)*strand
pick<-disttss<=1500&disttss>=-1500
#pick<-disttss<1e4&disttss>-1e4
disttssf<-disttss[pick]
logpf<-logp[pick]
############## Negative correlations###############
pick<-eqtl$b1<0
nlogp<-(-log10(as.numeric(eqtl$b1_p[pick])))
nmpos<-as.numeric(eqtl$mPos[pick])
ntss<-as.numeric(eqtl$ePos[pick])
nstrand<-as.integer(eqtl$strand[pick])
ndisttss<-(nmpos-ntss)*nstrand
pick<-ndisttss<=1500&ndisttss>=-1500
ndisttssf<-ndisttss[pick] 
nlogpf<-nlogp[pick]
pdf("Plot_1500bp.pdf",width=10,height=8)
plot(ndisttssf,nlogpf,col="red",pch=20,cex=.6,
     ylim=c(min(c(logpf,nlogpf)),max(c(logpf,nlogpf))),
     xlab="Distance between CpG sites and TSS (bp)",
     ylab="-log(P-value)",
     xaxt="n"
)
axis(1,at=seq(-1e4,1e4,length.out=5),labels=c("-100 00","-50 00", "0", "50 00","100 00"))
legend("topleft", c("Red: Negative correlations","Blue: Positive correlations"), text.col=c("red","blue"),bty="n", col=c("red","blue"), lwd=1, lty=c(0,0),pch=c(20,20))
points(disttssf,logpf,col="blue",pch=20,cex=.6)
dev.off()
#############################################################
# Quantsmooth: eQTL significance level
annot<-eqtl[,c(10,11)] # chromosome and position
colnames(annot)<-c("CHR","MapInfo")
annot$CHR[annot$CHR==23]<-"X"
annot$CHR[annot$CHR==24]<-"Y"
pdf("DistributionPlot_1500bp.pdf",width=10,height=8)
out<-prepareGenomePlot(annot,sexChromosomes=F,organism="hsa",units="bases",paintCytobands=T)
out<-data.frame(out)
logp<--log(eqtl$b1_p)
x<-out$MapInfo[eqtl$b1<0]
y<-(out$CHR+(logp/100))[eqtl$b1<0]
points(x=x,y=y,pch=20,cex=.4,col="red")
x<-out$MapInfo[eqtl$b1>0]
y<-(out$CHR+(logp/100))[eqtl$b1>0]
points(x=x,y=y,pch=20,cex=.5,col="blue")
dev.off()
save.image("eMAP_1500bp_CHOL.RData")
