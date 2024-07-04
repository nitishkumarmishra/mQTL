x <- read.table(file ="CpG-Island-HG19.tsv", header = T)
##m <- as.matrix(read.table(file="Matrix.csv", header=TRUE, sep = ",", row.names = 1, as.is=TRUE)) (Read CSV file in Matrix for MatrixStats)
x1 <-subset(x, obsExp >0.65) ## This command will take only those lines where ObsExp is greater thsn0.65
x2 <-x1[order(x1$obsExp) , ] ## Sort X1 in ascending order (Low to high )
x2 <-x1[order(-x1$obsExp) , ] ## sort x1 in descending order (High to low)
CpGislandDensity <- makeGRangesFromDataFrame(x1, keep.extra.columns = TRUE) ## This command convert X1 in genomic ranges and keeping all colum of data frame X1 in Iranges
upstream.probes <- names(subsetByOverlaps(hm450, flank(CpGislandDensity, 3000)))
downstream.probes<-names(subsetByOverlaps(hm450,flank(CpGislandDensity,-3000,start=TRUE)))
Island.3KbUpDown.probes<-unique(c(upstream.probes,downstream.probes))