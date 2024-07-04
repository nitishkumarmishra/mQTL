###################################################################
#### genomewide methylation profiling 
## This code is based on Figure 3 of Benton et al. Genome Biology (2015) An analysis of DNA methylation in human adipose tissue reveals differential modification of obesity genes before and after gastric bypass and weight loss

# I can't remember if all of the below packages are required, but I included them just in case
require(ggplot2)
require(ggbio)
require(biomaRt)
require(GenomeGraphs)
require(GenomicRanges)
library(BSgenome)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg19)

# get chromosome information
chr.len = seqlengths(Hsapiens)  # get chromosome lengths
chrom.length = chr.len[grep("_|M|Y", names(chr.len), invert = T)] # remove X,Y,M and random chromosomes if required

# create the ideogram accordingly
myIdeo <- GRanges(seqnames = names(chrom.length), ranges = IRanges(start = 1, width = chrom.length))
seqlevels(myIdeo) = names(chrom.length)
seqlengths(myIdeo) = (chrom.length)

# create object with genomic locations of methylation sites, as well as required statistics
map_data <- ab_data[c(1,5,5,4,11,12,13)] # subset data to include everything required
colnames(map_data)[2] <- "start"
colnames(map_data)[3] <- "end"
map_data$CHR <- paste("chr", map_data$CHR, sep="") # ensure chromosome labels are correct

# generate the 'inner' and 'outer' objects to plot (hypo and hyper methlation)
increase <- map_data[map_data$median_diff >= 0,]
decrease <- map_data[map_data$median_diff <= 0,]
increase <- map_data[map_data$median_diff >= 0,]
decrease <- map_data[map_data$median_diff <= 0,]
increase$strand = "*"
decrease$strand = "*"
colnames(increase)[4] <- "chr"
colnames(decrease)[4] <- "chr"
increase <- na.omit(increase) # ensure no NA values
decrease <- na.omit(decrease) # ensure no NA values
# can set a genome wide threshold if you don't want to plot every CpG site
decrease2 <- decrease[decrease$t.pval <= 1e-7,]
increase2 <- increase[increase$t.pval <= 1e-7,]

# create the granges objects
g.per.inc <- with(increase2, GRanges(chr, IRanges(start, end), strand = strand))
g.per.inc = keepSeqlevels(g.per.inc, names(chrom.length))
values(g.per.inc)$id = "increase"
values(g.per.inc)$p.val = increase2$t.pval
values(g.per.inc)$t.stat = increase2$t.stat
g.per.inc@seqinfo@seqlengths <- seqlengths(myIdeo)
#
g.per.dec <- with(decrease2, GRanges(chr, IRanges(start, end), strand = strand))
g.per.dec = keepSeqlevels(g.per.dec, names(chrom.length))
values(g.per.dec)$id = "decrease"
values(g.per.dec)$p.val = decrease2$t.pval
values(g.per.dec)$t.stat = decrease2$t.stat
g.per.dec@seqinfo@seqlengths <- seqlengths(myIdeo)

# circle plot using ggplot
p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", 
                              radius = 36, trackWidth = 3)

p <- p + layout_circle(c(g.per.inc, g.per.dec), geom = "point", size = 2.35,  aes(x = midpoint, y = t.stat, color = id), 
                       radius = 24, trackWidth = 30) + scale_colour_manual(values = c("deepskyblue4", "darkred")) 

p$labels$colour = "Differential Methylation" # label key accordingly

p + layout_circle(myIdeo, geom = "text", size = 5, face = "bold", aes(label = seqnames), 
                  vjust = 0, radius = 52, trackWidth = 6) + 
  theme(title = element_text("ESR Ab methylation analysis"), 
        legend.text = element_text(colour="black", size = 14),
        legend.title = element_text(colour="black", size=14)) 
#### END