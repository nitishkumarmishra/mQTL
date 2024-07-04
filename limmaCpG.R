########### limma for differential CpG sites Identification ##############
### Not CpG island #########
y <- c(rep(2,184),rep(1,10))
#design1 <- cbind(Grp1=2,TumorvsNormal=y, row.names=colnames(t1))
design1 <- data.frame(cbind(Grp1=2,TumorvsNormal=y), row.names=colnames(t1))## This is direct method to make contrast matrix
DMRfit <- lmFit(t1, design1)
DMRfitEb <- eBayes(DMRfit)
cutoff <- 0.01
DMR <- topTable(DMRfitEb, coef = "TumorvsNormal", number = Inf, p.value = cutoff, adjust.method="BH", sort.by = "p")
# head(DMR)
# logFC   AveExpr         t      P.Value    adj.P.Val        B
# cg08858649  0.3701431 0.4511610  8.531386 2.817804e-15 3.660610e-11 24.06908
# cg15834072  0.3073478 0.4557357  7.741638 4.021670e-13 2.228013e-09 19.31497
# cg01294808  0.3074352 0.5551836  7.659289 6.650973e-13 2.228013e-09 18.83325
# cg05336395  0.3073583 0.5261691  7.654204 6.860173e-13 2.228013e-09 18.80360
# cg14356919 -0.3235928 0.3377828 -7.579719 1.078454e-12 2.802040e-09 18.37049
# cg17902007  0.1478756 0.9623049  7.537850 1.389258e-12 3.007975e-09 18.12808
######################
region <- c(rep("Tumor", 184), rep("Normal", 10))
REG <- factor(region)
design <- model.matrix(~0+REG)
colnames(design)<-levels(REG)
fit <- lmFit(t1, design)
cont.matrix <- makeContrasts(TumorVsNormal=Tumor-Normal, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
cutoff <- 0.01
limres <- topTable(fit2, coef = 1, number = INF, p.value = cutoff, adjust.method = "BH", sort.by = "p")
# head(limres)
# logFC   AveExpr         t      P.Value    adj.P.Val        B
# cg08858649  0.3701431 0.4511610  8.531386 2.817804e-15 3.660610e-11 24.06908
# cg15834072  0.3073478 0.4557357  7.741638 4.021670e-13 2.228013e-09 19.31497
# cg01294808  0.3074352 0.5551836  7.659289 6.650973e-13 2.228013e-09 18.83325
# cg05336395  0.3073583 0.5261691  7.654204 6.860173e-13 2.228013e-09 18.80360
# cg14356919 -0.3235928 0.3377828 -7.579719 1.078454e-12 2.802040e-09 18.37049
# cg17902007  0.1478756 0.9623049  7.537850 1.389258e-12 3.007975e-09 18.12808
mt.rawp2adjp(DMR$P.Value, proc=c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY","ABH","TSBH"), alpha=0.01)
#############################################
################ CpG-Island analysis ########################
CpGcginame <- rownames(cginame)
matrixRowname <- row.names(matrix)
intersectProbe <- intersect(CpGcginame, matrixRowname)
beta.inCGI <- matrix[intersectProbe, ]
beta.CGI <- aggregate(beta.inCGI, by = list(intersectProbe), mean, na.rm = T)