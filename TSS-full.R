##############Biomart convert HGNC_SYMBOL IN ENTREZ-GENE-ID##########
library(biomaRt)
# define biomart object
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# read in the file
genes <- read.csv("TSS-HGNC-name.csv") ### TSS-HGNC.CSV is the list of TCGA genes in RNASeqV2 level-3 data
# query biomart
results <- getBM(attributes = c("entrezgene", "hgnc_symbol"), filters = "hgnc_symbol", values = genes$hgnc, mart = mart)
write.csv(results, file="biomart.csv")
x<-read.csv("biomart.csv")
yy<-as.character(x$entrezgene, na.rm=FALSE)
################ Function to call ENTREZ-GENE-ID and get 1.5kb up and downstream probes from TSS ########
library(Homo.sapiens)
keytypes(Homo.sapiens)
txs <- transcriptsBy(Homo.sapiens, 'gene', col='GENEID')
library(FDb.InfiniumMethylation.hg19)
hm450 <- get450k()
#use Entrez GeneID as a string, not numeric or factor

############################ Final line to get list of probes for all genes ##########
listprobes <- unlist(sapply(yy, getProbes))
#listprobes <- sapply(y, getProbes)
listprobes1 <-lapply(y, getProbes)
