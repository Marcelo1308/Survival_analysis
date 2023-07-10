
setwd("C:/Users/shiva/Desktop/Survival-workshop")

#Download gene expression data from TCGA using TCGA-Assembler2
geneExp2 <- DownloadRNASeqData(cancerType = "BRCA", assayPlatform
                               = "gene.normalized_RNAseq",tissueType = "TP",
                               saveFolderName = ".")
#Download clinical data from TCGA

setwd("C:/Users/shiva/Desktop/Survival-workshop")
getwd()
TCGAGEdata <- read.csv(file="BRCA__gene.normalized_RNAseq__TP__20230211101034.txt",sep = "\t",row.names = 1, check.names = F)
#TCGAMergedataOriginal <- read.csv(file ="TCGA-2.csv", sep = "\t",check.names = F,row.names = 1)
#dim(TCGAMergedataOriginal)

TCGAMergedata <- read.csv(file ="TCGA-3.csv", sep = "\t",check.names = F,row.names = 1)
dim(TCGAMergedata)
View(TCGAMergedata)
getwd()
table(TCGAMergedata$vital_status)
table(TCGAMergedata$`Overall survival`)
table(TCGAMergedata$race)
#surv_object <- Surv(time = TCGAMergedata$`Overall survival`, event = TCGAMergedata$vital_status)
#fit <- survfit(surv_object~TCGAMergedata$gender)

######Estimation of survival difference between white and black or African american
dataRace=subset(TCGAMergedata, (TCGAMergedata$race!='not reported')&(TCGAMergedata$race!='american indian or alaska native')&(TCGAMergedata$race!='asian'))
table(dataRace$race)
surv_object <- Surv(time = dataRace$`Overall survival`, event = dataRace$vital_status)
fit <- survfit(surv_object~dataRace$race) 
ggsurvplot(fit, data=dataRace,pval=TRUE)
surv_diff1 <- survdiff(Surv(dataRace$`Overall survival`, dataRace$vital_status) ~ dataRace$race, data = dataRace)
surv_diff1
ggsurvplot(fit,data= dataRace,
           pval=TRUE, conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

ggsurvplot(fit,data=dataRace,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##Estimation of survival difference between early and late stages of breast cancer
table(TCGAMergedata$ajcc_pathologic_stage)
dataStage=mutate(TCGAMergedata, ajcc_pathologic_stage=ifelse(((ajcc_pathologic_stage=="Stage I")|(ajcc_pathologic_stage=="Stage IA")|(ajcc_pathologic_stage=="Stage IB")|(ajcc_pathologic_stage=="Stage II")|(ajcc_pathologic_stage=="Stage IIA")|(ajcc_pathologic_stage=="Stage IIA")|(ajcc_pathologic_stage=="Stage IIB")),"Early-Stage","Late-Stage"))
table(dataStage$ajcc_pathologic_stage)
surv_object1 <- Surv(time = dataStage$`Overall survival`, event = dataStage$vital_status)
fit <- survfit(surv_object1~dataStage$ajcc_pathologic_stage) 
ggsurvplot(fit, data=dataStage,pval=TRUE)
surv_diff2 <- survdiff(Surv(dataStage$`Overall survival`, dataStage$vital_status) ~ dataStage$ajcc_pathologic_stage, data = dataStage)
surv_diff2
ggsurvplot(fit,data= dataStage,
           pval=TRUE, conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

ggsurvplot(fit,data=dataStage,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

####Effect of gene expression in survival#BRAF
surv_object2 <- Surv(time = TCGAMergedata$`Overall survival`, event = TCGAMergedata$vital_status)
which(is.na(TCGAMergedata$`BRAF|673`))
fit3 <- survfit(surv_object2~TCGAMergedata$`BRAF|673`>median(TCGAMergedata$`BRAF|673`))
ggsurvplot(fit3, data=TCGAMergedata,pval=TRUE)
#surv_diff3 <- survdiff(Surv(TCGAMergedata$`Overall survival`, TCGAMergedata$vital_status) ~ TCGAMergedata$`BRAF|673`, data = TCGAMergedata)
#surv_diff3

#Effect of gene expression on survival##BRCA1
which(is.na(TCGAMergedata$`TP53|7157`))
fit5 <- survfit(surv_object2~TCGAMergedata$`BRCA1|672`>median(TCGAMergedata$`BRCA1|672`))
ggsurvplot(fit5, data=TCGAMergedata,pval=TRUE)
