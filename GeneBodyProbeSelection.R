## R code for getting CpG's from Gene body
## In this code we also run some linux commands
library(Homo.sapiens)
#keytypes(Homo.sapiens)
txs <- transcriptsBy(Homo.sapiens, 'gene', col='GENEID')
library(FDb.InfiniumMethylation.hg19)

hm450 <- get450k()

getProbes<-function(geneID){
  temp<-txs[[geneID]]
  if(is.null(temp)){
    probes<-NULL
  }else{
    upstream.probes<-names(subsetByOverlaps(hm450,temp))
    probes<-unique(c(upstream.probes))
  }
  return(probes)
}
GeneID <- names(txs)
yy <- GeneID[1:10]
listprobes <- sapply(yy, getProbes)



tmp <- listprobes[lapply(listprobes,length)>0] ## Remove character(0) from list
tmp <- tmp[!sapply(tmp, is.null)] ## Remove NULL from list
names <- names(tmp)# may be its not necessary, because tmp can give rownames
tmp <- do.call(rbind, listprobes)
write.csv(tmp, file="tmp.csv")

############
dim(tmp)# 1292 column 
########### Linux command to make file
sed -ie 's/","/\t/g' tmp.csv
sed -ie 's/"//g' tmp.csv
cut -d ' ' -f 1 tmp.csv > tt1
cut -d ' ' -f 2-1292 tmp.csv > tt2
sed -ie 's/      /;/g' tt2
paste -d "       " tt1 tt2  >tmp1.csv
awk '{ split($2,a,";"); for (i in a) print $1 "\t" a[i]; }' tmp1.csv|sort -u >ProbesGeneBody.txt
##############
GeneBody <- read.csv("ProbesGeneBody.txt", header=FALSE, sep="\t")
GeneBody$V1 <- as.character(GeneBody$V1)
GeneBody$V2 <- as.character(GeneBody$V2)
colnames(GeneBody) <- c("GeneID", "ProbeID")

save.image("ProbesGeneBody.RData")
