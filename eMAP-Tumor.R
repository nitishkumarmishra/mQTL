library(eMap) # see references in paper for web page
library(quantsmooth)
# --- Illumina annotation
# http://support.illumina.com/downloads/humanmethylation450_15017482_v1-2_product_files.html
#data<-read.csv("HumanMethylation450_15017482_v1-2.csv",header=T,sep=",",na.strings="",row.names=NULL,quote="",skip=7)
#data<-read.csv("HumanMethylation450_15017482_v1-2.csv",header=T,sep=",",na.strings="",row.names=NULL,quote="",skip=7)
#data<-data[-(485578:nrow(data)),] #remove controls
#probeinfo<-data[,colnames(data)%in%c("IlmnID","CHR","MAPINFO","Relation_to_UCSC_CpG_Island","Enhancer","Infinium_Design_Type","UCSC_CpG_Islands_Name")]
#colnames(probeinfo)<-c("Probe","Design","Chr","Pos","CGI_name","Relation_CGI","Enhancer")
#write.table(probeinfo, file="Probeinfo.txt", row.names=F, quote=F, sep="\t")


# --- Make table of genes containing gene name, TSS and strand. 

# Go to UCSC table browser (http://genome.ucsc.edu/cgi-bin/hgTables)
# assembly: hg19
# group: Genes and Gene Predictions
# track: RefSeq Genes
# table: refGenes
# output format: all fields from selected table
# output file: refseq_genes_hg19.txt

#data<-read.table("refseq_genes_hg19.txt",header=T,sep="\t",na.strings="",row.names=NULL,quote="",comment.char="")
#data<-data[,colnames(data)%in%c("name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","name2")]

# The following code removes alternative transcripts for each gene, so that each gene only has one TSS. 
#res<-data[1,]
#k<-0
#for(g in unique(data$name2)){
#	k<-k+1;print(k)
#	d<-data[data$name2==g,]
#	if(nrow(d)==1){res<-rbind(res,d)}
#	if(nrow(d)>1){
#		d<-d[grep("NM",d$name),]
#		d$char<-as.numeric(lapply(as.character(d$name),nchar))
#		d$num<-as.numeric(gsub("NM_","",d$name))
#		d<-d[order(d$char,d$num),]
#		d<-d[1,1:ncol(res)] # This approach will chose the NM number with 1) the fewest characters, and, if equal 2) the lowest number. This might not be ideal, but it't the best I could think of.
#		res<-rbind(res,d)
#	}
#}
#res<-res[-1,]
#res<-res[!is.na(res$name2),]
#write.table(res, file="refseq_genes_hg19_finished.txt", row.names=F, quote=F, sep="\t")

# Now you know the TSS of all genes, and from Illumina annotation you know the location of the probes


# --- eMap -- meth-expression analysis
#expr<-read.table("...",header=T,sep="\t",na.strings="NA",row.names=NULL)
	# This file should contain 4 information columns with "Gene", "Chr", "TSS" and "Strand". The rest of the columns should be expression data of samples. 
	# or make the object in R.
res <- read.csv("refseq_genes_hg19_finished.txt", header=TRUE, sep="\t")
probeinfo <- read.csv("Probeinfo.txt", header=TRUE, sep="\t")
load("GeneExp.Rda")
load("BMIQ_Meth.Rda")
expr <- round(geneExp)
expr <- expr[rowSums(expr)>0,]
expr <- log2(expr+1)
expr <- expr[,grep("01A", colnames(expr))]
expr$gene <- sapply(strsplit(rownames(expr),"\\|"),'[[',1)
exprMerge <- merge(expr, res, by.x="gene", by.y="name2")
exprMerge$chrom <- gsub("chr","", exprMerge$chrom)
#exprMerge$chrom <- gsub("chr","", exprMerge$chrom)
colnames(exprMerge)[39] <- c("Chr")
colnames(exprMerge)[41] <- c("TSS")
expr <- exprMerge[order(exprMerge$Chr,exprMerge$TSS),]
expr<-expr[order(expr$Chr,expr$TSS),]
#expr <- expr[-which(expr$Chr=="Un_gl000228"),]
me<-as.matrix(expr[,2:37]);rownames(me)<-expr[,1]
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
#me <- me[,grep("01A", colnames(me))]
#mm <- mm[,grep("01A", colnames(mm))]
eMap1(me=me,mm1=mm,output.tag="output",p.cut=1,cis.only=T,cis.distance=1e5,eChr=eChr,ePos=ePos,mChr=mChr,mPos=mPos)
eqtl<-read.table("output_eqtl.txt",header=T,sep="\t",na.strings="NA",row.names=NULL)
eqtl<-eqtl[,colnames(eqtl)%in%c("Gene_ID","Marker_ID","b1","b1_p")]
eqtl$bonf_p<-p.adjust(eqtl$b1_p,method="bonferroni")
eannot<-cbind(data.frame(nr=1:nrow(me)),Gene=rownames(me),eChr=eChr,ePos=ePos)
mannot<-cbind(data.frame(nr=1:nrow(mm)),Probe=rownames(mm),mChr=mChr,mPos=mPos)
eqtl<-merge(eqtl,eannot,all.x=T,all.y=F,by.x="Gene_ID",by.y="nr")
eqtl<-merge(eqtl,mannot,all.x=T,all.y=F,by.x="Marker_ID",by.y="nr")

#eqtl<-merge(eqtl,expr[,colnames(expr)%in%c("Gene","Strand")],all.x=T,all.y=F,by.x="Gene",by.y="Gene")
eqtl <- merge(eqtl,expr[,colnames(expr)%in%c("gene","strand")],all.x=T,all.y=F,by.x="Gene",by.y="gene")
#eqtl<-merge(eqtl,me2[,colnames(me2)%in%c("Gene","Strand")],all.x=T,all.y=F,by.x="Gene",by.y="Gene")
#eqtl$strand <- gsub("\\+", "+1", eqtl$strand); eqtl$strand <- gsub("\\-", "+1", eqtl$strand)
eqtl$dist<-(eqtl$mPos-eqtl$ePos)*eqtl$strand
saveRDS(eqtl,"output_eqtl_annot.rds")
#eqtl<-eqtl[eqtl$bonf_p<0.05,]
eqtl <- eqtl[which(eqtl$bonf_p<0.05),]
saveRDS(eqtl,"output_eqtl_annot_bonf.rds")
write.csv(eqtl, file="eqtl-0.5.csv")
# Plot
# Positive correlations
# pick<-eqtl$b1>0
# logp<-(-log10(as.numeric(eqtl$b1_p[pick])))
# mpos<-as.numeric(eqtl$mPos[pick])
# tss<-as.numeric(eqtl$TSS[pick])
# strand<-as.integer(eqtl$Strand[pick])
# disttss<-(mpos-tss)*strand
# pick<-disttss<1e5&disttss>-1e5
# disttssf<-disttss[pick]
# logpf<-logp[pick]
# # Negative correlations
# pick<-eqtl$b1<0
# nlogp<-(-log10(as.numeric(eqtl$b1_p[pick])))
# nmpos<-as.numeric(eqtl$mPos[pick])
# ntss<-as.numeric(eqtl$TSS[pick])
# nstrand<-as.integer(eqtl$Strand[pick])
# ndisttss<-(nmpos-ntss)*nstrand
# pick<-ndisttss<1e5&ndisttss>-1e5
# ndisttssf<-ndisttss[pick] 
# nlogpf<-nlogp[pick]
# 
# pdf("plot.pdf",width=10,height=8)
# plot(ndisttssf,nlogpf,col="red",pch=20,cex=.6,
# 	ylim=c(min(c(logpf,nlogpf)),max(c(logpf,nlogpf))),
# 	xlab="Distance between methylation probe and TSS (bp)",
# 	ylab="-log p-value",
# 	xaxt="n"
# 	)
# axis(1,at=seq(-1e5,1e5,length.out=5),labels=c("-100 000","-50 000","0","50 000","100 000"))
# legend("topleft", c("Red: negative correlations","Blue: positive correlations"), text.col=c("red","blue"),bty="n")
# points(disttssf,logpf,col="blue",pch=20,cex=.6)
# dev.off()
#############################################################
###### Positive correlations ############
pick<-eqtl$b1>0
logp<-(-log10(as.numeric(eqtl$b1_p[pick])))
mpos<-as.numeric(eqtl$mPos[pick])
tss<-as.numeric(eqtl$ePos[pick])
strand<-as.integer(eqtl$strand[pick])
disttss<-(mpos-tss)*strand
pick<-disttss<1e5&disttss>-1e5
#pick<-disttss<1e4&disttss>-1e4
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

pdf("Plot_cancer.pdf",width=10,height=8)
plot(ndisttssf,nlogpf,col="red",pch=20,cex=.6,
     ylim=c(min(c(logpf,nlogpf)),max(c(logpf,nlogpf))),
     xlab="Distance between CpG sites and TSS (bp)",
     ylab="-log(P-value)",
     xaxt="n"
)

#axis(1, xaxp=c(-100000, 100000, 5), las=2)
axis(1,at=seq(-1e5,1e5,length.out=5),labels=c("-50 000","-25 000", "0", "25 000","50 000"))
#axis(1,at=seq(-1e4,1e4,length.out=5),labels=c("-10 000","-5 000","0","5 000","10 000"))
legend("topleft", c("Red: Negative correlations","Blue: Positive correlations"), text.col=c("red","blue"),bty="n", col=c("red","blue"), lwd=1, lty=c(0,0),pch=c(20,20))
points(disttssf,logpf,col="blue",pch=20,cex=.6)
dev.off()
save.image("CHOLemap_cancer.RData")

#############################################################
# --- Quantsmooth: eQTL significance level
annot<-eqtl[,c(10,11)] # chromosome and position
colnames(annot)<-c("CHR","MapInfo")
annot$CHR[annot$CHR==23]<-"X"
annot$CHR[annot$CHR==24]<-"Y"
# annot <- annot[,-(annot$CHR==23)]
# annot <- annot[,-(annot$CHR==24)]
#bmp("ChromosomeDistribution.jpg",height=1200,width=1200)
#bmp("ChromosomeDistribution.jpg",height=500,width=500)
pdf("Plot.pdf",width=10,height=8)
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








