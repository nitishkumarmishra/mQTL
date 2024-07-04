## Make clinical and GeneExp data in same order
GeneExp <- GeneExp[,order(substr(colnames(GeneExp), 1,12))]
Clinic1$status <- ifelse(substr(Clinic1$sampleID, 14,15)=="01" , "Normal", "Tumor")
Clinic1$cancer_types <- ifelse(Clinic1$histologic_diagnosis=="Cholangiocarcinoma; intrahepatic" , "intrahepatic", "extrahepatic")
design <- model.matrix(~0+satus+types, data = Clinic1)
colnames(design) <- c(levels(satus), levels(types)[-1])
