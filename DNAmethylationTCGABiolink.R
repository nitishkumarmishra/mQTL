library(TCGAbiolinks)
library(SummarizedExperiment)

###################
library("preprocessCore")
normalize.quantiles(datamatrix,copy=TRUE)
## do quantile normalizatio of TCGA methylation data 
## although level3 is normalized but do quantile normalization

#Meth is TCGA level3 beta value file
champ.norm(beta = Meth, resultsDir = paste(getwd(), "resultsChamp",sep = "/"), methValue = "B", fromIDAT = FALSE, norm = "BMIQ", fromFile = TRUE, betaFile, filter = FALSE, filterXY = TRUE, QCimages = FALSE, plotBMIQ = FALSE)


######## DNA methylation ################
cancer <- "TGCT"
dataType <- "HumanMethylation450"
pathCancer <- paste0("data",cancer)
PlatformCancer <- "HumanMethylation450"
datQuery <- TCGAquery(tumor = cancer, platform = PlatformCancer, level = "3")
lsSample <- TCGAquery_samplesfilter(query = datQuery)
dataSmTP <- TCGAquery_SampleTypes(barcode = lsSample$HumanMethylation450,
                                  typesample = "TP")
dataSmTN <- TCGAquery_SampleTypes(barcode = lsSample$HumanMethylation450,
                                  typesample ="NT")
dataClin <- TCGAquery_clinic(tumor = cancer,
                             clinical_data_type = "clinical_patient")


## Get TCGA ID's 
which(dataClin$histological_type=="Seminoma; NOS")
paad.met <- TCGAprepare(query = datQuery,
                        dir = "dataTGCT",
                        save = TRUE,
                        filename = "metpaad.rda",
                        reannotate = FALSE)

getTCGA("TGCT", Meth = TRUE, RNA = FALSE, Clinic = FALSE, basedir = "./dataTGCT/", RNAtype = "gene", Methfilter = 0.2)
## getTCGA from ELMER
