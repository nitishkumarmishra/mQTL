library(eMap) # see references in paper for web page
library(quantsmooth)
res <- read.csv("refseq_genes_hg19_finished.txt", header=TRUE, sep="\t")
probeinfo <- read.csv("Probeinfo.txt", header=TRUE, sep="\t")
load("GeneExp.Rda")
load("BMIQ_Meth.Rda")
expr <- round(geneExp)
expr <- expr[rowSums(expr)>0,]
expr <- log(expr+1)
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
save.image("CHOLeMAP.RData")
### eMAP analysis ################
eMap1(me=me,mm1=mm,output.tag="output",p.cut=1,cis.only=T,cis.distance=1e5,eChr=eChr,ePos=ePos,mChr=mChr,mPos=mPos)
eqtl<-read.table("output_eqtl.txt",header=T,sep="\t",na.strings="NA",row.names=NULL)
eqtl<-eqtl[,colnames(eqtl)%in%c("Gene_ID","Marker_ID","b1","b1_p")]
eqtl$bonf_p<-p.adjust(eqtl$b1_p,method="bonferroni")
eannot<-cbind(data.frame(nr=1:nrow(me)),Gene=rownames(me),eChr=eChr,ePos=ePos)
mannot<-cbind(data.frame(nr=1:nrow(mm)),Probe=rownames(mm),mChr=mChr,mPos=mPos)
eqtl<-merge(eqtl,eannot,all.x=T,all.y=F,by.x="Gene_ID",by.y="nr")
eqtl<-merge(eqtl,mannot,all.x=T,all.y=F,by.x="Marker_ID",by.y="nr")
eqtl <- merge(eqtl,expr[,colnames(expr)%in%c("gene","strand")],all.x=T,all.y=F,by.x="Gene",by.y="gene")
eqtl$dist<-(eqtl$mPos-eqtl$ePos)*eqtl$strand
saveRDS(eqtl,"output_eqtl_annot.rds")
eqtl<-eqtl[eqtl$bonf_p<0.05,]
saveRDS(eqtl,"output_eqtl_annot_bonf.rds")
write.csv(eqtl, file="eqtl-0.5.csv")
############# Positive correlations ##############
pick<-eqtl$b1>0
logp<-(-log10(as.numeric(eqtl$b1_p[pick])))
mpos<-as.numeric(eqtl$mPos[pick])
tss<-as.numeric(eqtl$ePos[pick])
strand<-as.integer(eqtl$strand[pick])
disttss<-(mpos-tss)*strand
pick<-disttss<1e5&disttss>-1e5
disttssf<-disttss[pick]
logpf<-logp[pick]
####### Negative correlations############
pick<-eqtl$b1<0
nlogp<-(-log10(as.numeric(eqtl$b1_p[pick])))
nmpos<-as.numeric(eqtl$mPos[pick])
ntss<-as.numeric(eqtl$ePos[pick])
nstrand<-as.integer(eqtl$strand[pick])
ndisttss<-(nmpos-ntss)*nstrand
pick<-ndisttss<1e5&ndisttss>-1e5
ndisttssf<-ndisttss[pick] 
nlogpf<-nlogp[pick]

pdf("Plot.pdf",width=10,height=8)
plot(ndisttssf,nlogpf,col="red",pch=20,cex=.6,
     ylim=c(min(c(logpf,nlogpf)),max(c(logpf,nlogpf))),
     xlab="Distance between CpG sites and TSS (bp)",
     ylab="-log(P-value)",
     xaxt="n"
)

#axis(1, xaxp=c(-100000, 100000, 5), las=2)
axis(1,at=seq(-1e5,1e5,length.out=5),labels=c("-100 000","-50 000", "0", "50 000","100 000"))
#axis(1,at=seq(-1e4,1e4,length.out=5),labels=c("-10 000","-5 000","0","5 000","10 000"))
legend("topleft", c("Red: Negative correlations","Blue: Positive correlations"), text.col=c("red","blue"),bty="n", col=c("red","blue"), lwd=1, lty=c(0,0),pch=c(20,20))
points(disttssf,logpf,col="blue",pch=20,cex=.6)
dev.off()


#############################################################
# --- Quantsmooth: eQTL significance level
annot<-eqtl[,c(10,11)] # chromosome and position
colnames(annot)<-c("CHR","MapInfo")
annot$CHR[annot$CHR==23]<-"X"
annot$CHR[annot$CHR==24]<-"Y"
pdf("DistributionPlot.pdf",width=10,height=8)
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

save.image("eMAP_logNormRead.RData")






