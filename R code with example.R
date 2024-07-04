
library(Homo.sapiens)

keytypes(Homo.sapiens)

txs <- transcriptsBy(Homo.sapiens, 'gene', col='GENEID')


library(FDb.InfiniumMethylation.hg19)

hm450 <- get450k()

#use Entrez GeneID as a string, not numeric or factor
                                        
getProbes<-function(geneID){
  temp<-txs[[geneID]]
if(is.null(temp))
{

  probes<-NULL

}
else
{
upstream.probes<-names(subsetByOverlaps(hm450,flank(temp,1500)))
downstream.probes<-names(subsetByOverlaps(hm450,flank(temp,-1500,start=TRUE)))
probes<-unique(c(upstream.probes,downstream.probes))

}
return(probes)
}
#example
getProbes('1234')
