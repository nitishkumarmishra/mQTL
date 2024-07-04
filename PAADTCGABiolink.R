library(TCGAbiolinks)
library(SummarizedExperiment)
######## DNA methylation ################
cancer <- "PAAD"
dataType <- "HumanMethylation450"
pathCancer <- paste0("data",cancer)

query.met <- TCGAquery(tumor = c("paad"),
                       platform = c("HumanMethylation450"), level = 3)
TCGAdownload(query.met, path = pathCancer)
paad.met <- TCGAprepare(query = query.met,
                        dir = "dataPAAD",
                        save = TRUE,
                        filename = "metpaad.rda",
                        reannotate = TRUE)
paad.met <- subset(paad.met,select = !(colData(paad.met)$methylation_subtype %in% c(NA)))
TCGAvisualize_meanMethylation(paad.met,
                              groupCol = "methylation_subtype",
                              subgroupCol = "hypermethylated",
                              group.legend  = "Groups",
                              subgroup.legend = "hypomethylated",
                              filename = "paad_mean.png")

########### Gene expression #############
cancer <- "PAAD"
dataType <- "rsem.genes.results"
pathCancer <- paste0("data",cancer)
PlatformCancer <- "IlluminaHiSeq_RNASeqV2"
datQuery <- TCGAquery(tumor = cancer, platform = PlatformCancer, level = "3")
lsSample <- TCGAquery_samplesfilter(query = datQuery)
dataSmTP <- TCGAquery_SampleTypes(barcode = lsSample$IlluminaHiSeq_RNASeqV2,
                                  typesample = "TP")
dataSmTN <- TCGAquery_SampleTypes(barcode = lsSample$IlluminaHiSeq_RNASeqV2,
                                  typesample ="NT")
dataClin <- TCGAquery_clinic(tumor = cancer,
                             clinical_data_type = "clinical_patient") 
TCGAdownload(datQuery, path = pathCancer, type = dataType, samples = c(dataSmTP, dataSmTN))
PAADnaseq_assay <- TCGAprepare(datQuery,"dataPAAD",type = "rsem.genes.results")
PAADMatrix <- assay(PAADnaseq_assay,"raw_counts")
PAADRnaseq_CorOutliers <- TCGAanalyze_Preprocessing(PAADnaseq_assay)

dataAssy <- TCGAprepare(query = datQuery,
                        dir = pathCancer,
                        type = dataType,
                        save = TRUE,
                        summarizedExperiment = TRUE,
                        samples = c(dataSmTP,dataSmTN),
                        filename = paste0(cancer,"_",PlatformCancer,".rda"))

dataPrep <- TCGAanalyze_Preprocessing(object = dataAssy, cor.cut = 0.6)  
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent") 
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut =  0.25)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmTN],
                            mat2 = dataFilt[,dataSmTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")  
ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",
                                RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = rownames(dataDEGs),
                        nBar = 20)
TCGAVisualize_volcano(dataDEGs$logFC,dataDEGs$FDR,
                      filename = "PAAD_volcanoexp.png",
                      x.cut = 2,
                      y.cut = 0.01,
                      names = rownames(dataDEGs),
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      title = "Volcano plot (Tumor vs Normal)")


#dataClin$neoplasm_histologic_grade
#dataClin$pathologic_T
#dataClin$pathologic_N
#dataClin$pathologic_stage
#dataClin$person_neoplasm_cancer_status