data<-read.table("TCGA_GDC_Harmonided_RNASeq-GFT.txt",header=T,sep="\t",na.strings="",row.names=NULL,quote="",comment.char="")
# http://www.uniprot.org/docs/pkinfam # Total 504 Human Kinase genes on August 30 2017
data1<-read.csv("Uniport_Kinase.txt",header=T,sep="\t",na.strings="",row.names=NULL) ## Total 504 Uniprot Kinase genes  
data1[data1$Symbol %in% data$Gene_Symbol,]$Symbol ### Total 496 Human Protein kinases in TCGA Harmonized data
data1[!data1$Symbol %in% data$Gene_Symbol,]$Symbol## 8 Kinase genes are not present in TCGA Harmonized


Kinase_absent <- data1[!data1$Symbol %in% data$Gene_Symbol,]$Symbol
Kinase_present <- data1[data1$Symbol %in% data$Gene_Symbol,]$Symbol
save.image(file = "Kinase_Uniprot.RData")
