
library(caOmicsV)

data(biomatrixPlotDemoData)

plotBioMatrix(biomatrixPlotDemoData, summaryType="text")

bioMatrixLegend(heatmapNames=c("RNASeq", "miRNASeq"),
                
                categoryNames=c("Methyl H", "Methyl L"),
                
                binaryNames=c("CN LOSS", "CN Gain"),
                
                heatmapMin= ???3, heatmapMax=3, colorType="BlueWhiteRed")




data(bionetPlotDemoData)

plotBioNetCircos(bionetPlotDemoData)

dataNames<- c("Tissue Type", "RNASeq", "miRNASeq", "Methylation", "CNV")

bioNetLegend(dataNames, heatmapMin= ???3, heatmapMax=3)
