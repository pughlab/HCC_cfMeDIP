# Description: This script was employed to analysis HCC subtypes based on the different liver disease backgrounds using “one-vs-each” tumor subtyping machine learning strategy.
## Author：Kui Chen → kui.chen@uhn.ca, Ping Luo → ping.luo@uhn.ca
## Date Created: 2022-09-01
## Last Modified: 2023-12-10
# Before starting, do the followings
# 1.using source() to introduce the step_7s_one_vs_each_hcc.R code
# 2.require "n236_lcpm.rds", "sample_124.rds"

## Abbreviations ####################################################################################
# ALD, alcoholic liver disease
# CTL, healthy control
# HBV, hepatitis B virus
# HCV, hepatitis C virus
# MASH, Metabolic dysfunction-associated steatohepatitis
# ND, no known liver disease ########################################################################

library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(caret)
library(randomForest)
library(glmnet)
library(limma)
library(doParallel)
Cluster <-makeCluster(20)
Cluster <- registerDoParallel(Cluster)

setwd("/path/to/your/counts")

# Prepare matrix
load("n236_lcpm.RData")

#Prepare matrix
count <- readRDS('n236_lcpm.rds')
phenotype <- readRDS('sample_124.rds')

table(phenotype$Disease)
table(phenotype$Cirr_group)

matrix <- count[ ,match(phenotype$Sample_ID,colnames(count))]
identical(colnames(matrix), phenotype$Sample_ID)

#check how many DMR in black_list
load("black_bin_v2.RData")
matrix <- subset(matrix, !rownames(matrix) %in% black_bin$black_bin)
identical(colnames(matrix), phenotype$Sample_ID)

source("step_7s_one_vs_each_hcc.R")

Features.CVparam <- trainControl(method = "repeatedcv",number = 10, repeats =3, verboseIter=TRUE,returnData=FALSE, classProbs = TRUE,savePredictions=FALSE)

# Define the subclasses
subclasses <- c("HBV", "HCV", "ALD", "MASH", "ND", "CTL")

# Create lists to store the results
All.kFold <- list()
TestPerformance.list <- list()
Features.list <- list()

for (subclass in subclasses) {
 
 Splits10 <- list()
 
 for (j in 1:100) {
  # Split the data for 100 times, result of each split is saved in Splits10
  Splits10[[j]] <- SplitkFold(matrix, phenotype$Disease, 10)
  
  # 10-fold cross validation to calculate the probability (score) of each sample being the subclass
  kFold.list <- list()
  for(i in 1:10) {
   kFold.list[[i]] <- OnevsEach.HCC(matrix, classes.df = Splits10[[j]]$df, Indices = Splits10[[j]]$samples[[i]], nDMR = 300, subclass)
  }
  
  Prediction.classProbs <- kFold.list[[1]]$TestPred
  for (i in 2:10) {
   Prediction.classProbs <- bind_rows(Prediction.classProbs, kFold.list[[i]]$TestPred)
  }
  
  All.kFold[[j]] <- kFold.list
  
  TestPerformance.list[[j]] <- Prediction.classProbs
  
 }
 
 saveRDS(All.kFold, paste0(subclass, '_kFoldModel.rds')) # This file saves all the results for the subclass
 saveRDS(TestPerformance.list, paste0(subclass, '_TestPerformance.rds')) # This file can be used by the original ROC codes for the subclass
 
 Features.onetime <- unique(unlist(lapply(All.kFold, function(x) unlist(lapply(x, function(y) y$Model$finalModel$xNames)))))
 Features.list[[subclass]] <- Features.onetime
 saveRDS(Features.list, paste0(subclass, '_Features.list.rds')) # This file can be used for DMR selection
 
 Features <- unique(unlist(Features.list))
 save(Features.list, All.kFold, TestPerformance.list,Features, file = paste0('hcc_', subclass, "_res.RData"))

}

#save(Features.list, All.kFold, TestPerformance.list, file = "hcc_subtype_res.RData") # The features for each subclass are saved here
#load("hcc_subtype_res.RData")

# Obtain AUROC
library(pROC)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(grid)

# Draw AUROC Curves
# HBV Subclass
load("hcc_HBV_res.RData")
hbv_predictions <- lapply(TestPerformance.list, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "HBV",1,0))%>%mutate(Class2 = as.numeric(Class2)))
for(i in 1:length(hbv_predictions)) { hbv_predictions[[i]]$Index <- paste0("Model_",i) }
ROCs.computed <- lapply(hbv_predictions, function(x) with(x, roc(Class2 ~ One)))
ROCs.computed <- lapply(ROCs.computed, function(x) smooth(x, method = "density"))
ROCs.computed <- lapply(ROCs.computed, function(x) data.frame(Sens = x$sensitivities, Spec = x$specificities, stringsAsFactors = F))
AUCs <- lapply(hbv_predictions, function(x) with(x, roc(Class2 ~ One)$auc))
AUCs <- unlist(AUCs)
# Obtain AUROC
t.test(AUCs, conf.level = 0.95)
# Combine all model predictions into a large data frame
combined_predictions <- do.call(rbind, hbv_predictions)
# Compute the overall ROC
roc_hbv <- with(combined_predictions, roc(Class2 ~ One))
# Compute the smoothed ROC
roc_hbv_smooth <- smooth(roc_hbv, method = "density")
# Extract the sensitivities and specificities
roc_hbv <- data.frame(Sens = roc_hbv_smooth$sensitivities,Spec = roc_hbv_smooth$specificities, stringsAsFactors = F)
# table(TestPerformance.list[[1]]$ActualClass)

# CTL Subclass
load("hcc_CTL_res.RData")
ctl_predictions <- lapply(TestPerformance.list, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "CTL",1,0))%>%mutate(Class2 = as.numeric(Class2)))
for(i in 1:length(ctl_predictions)) { ctl_predictions[[i]]$Index <- paste0("Model_",i) }
ROCs.computed <- lapply(ctl_predictions, function(x) with(x, roc(Class2 ~ One)))
ROCs.computed <- lapply(ROCs.computed, function(x) smooth(x, method = "density"))
ROCs.computed <- lapply(ROCs.computed, function(x) data.frame(Sens = x$sensitivities, Spec = x$specificities, stringsAsFactors = F))
AUCs <- lapply(ctl_predictions, function(x) with(x, roc(Class2 ~ One)$auc))
AUCs <- unlist(AUCs)
# Obtain AUROC
t.test(AUCs, conf.level = 0.95)
# Combine all model predictions into a large data frame
combined_predictions <- do.call(rbind, ctl_predictions)
# Compute the overall ROC
roc_ctl <- with(combined_predictions, roc(Class2 ~ One))
# Compute the smoothed ROC
roc_ctl_smooth <- smooth(roc_ctl, method = "density")

# Extract the sensitivities and specificities
roc_ctl <- data.frame(Sens = roc_ctl_smooth$sensitivities, Spec = roc_ctl_smooth$specificities, stringsAsFactors = F)

# ALD Subclass
load("hcc_ALD_res.RData")
ald_predictions <- lapply(TestPerformance.list, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "ALD",1,0))%>%mutate(Class2 = as.numeric(Class2)))
for(i in 1:length(ald_predictions)) { ald_predictions[[i]]$Index <- paste0("Model_",i) }
ROCs.computed <- lapply(ald_predictions, function(x) with(x, roc(Class2 ~ One)))
ROCs.computed <- lapply(ROCs.computed, function(x) smooth(x, method = "density"))
ROCs.computed <- lapply(ROCs.computed, function(x) data.frame(Sens = x$sensitivities, Spec = x$specificities, stringsAsFactors = F))
AUCs <- lapply(ald_predictions, function(x) with(x, roc(Class2 ~ One)$auc))
AUCs <- unlist(AUCs)
# Obtain AUROC
t.test(AUCs, conf.level = 0.95)
# Combine all model predictions into a large data frame
combined_predictions <- do.call(rbind, ald_predictions)
# Compute the overall ROC
roc_ald <- with(combined_predictions, roc(Class2 ~ One))
# Compute the smoothed ROC
roc_ald_smooth <- smooth(roc_ald, method = "density")
# Extract the sensitivities and specificities
roc_ald <- data.frame(Sens = roc_ald_smooth$sensitivities,  Spec = roc_ald_smooth$specificities, stringsAsFactors = F)


# MASH Subclass
load("hcc_MASH_res.RData")
MASH_predictions <- lapply(TestPerformance.list, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "MASH",1,0))%>%mutate(Class2 = as.numeric(Class2)))
for(i in 1:length(MASH_predictions)) { MASH_predictions[[i]]$Index <- paste0("Model_",i) }
ROCs.computed <- lapply(MASH_predictions, function(x) with(x, roc(Class2 ~ One)))
ROCs.computed <- lapply(ROCs.computed, function(x) smooth(x, method = "density"))
ROCs.computed <- lapply(ROCs.computed, function(x) data.frame(Sens = x$sensitivities, Spec = x$specificities, stringsAsFactors = F))
AUCs <- lapply(MASH_predictions, function(x) with(x, roc(Class2 ~ One)$auc))
AUCs <- unlist(AUCs)
# Obtain AUROC
t.test(AUCs, conf.level = 0.95)
# Combine all model predictions into a large data frame
combined_predictions <- do.call(rbind, MASH_predictions)
# Compute the overall ROC
roc_MASH <- with(combined_predictions, roc(Class2 ~ One))
# Compute the smoothed ROC
roc_MASH_smooth <- smooth(roc_MASH, method = "density")
# Extract the sensitivities and specificities
roc_MASH <- data.frame(Sens = roc_MASH_smooth$sensitivities,  Spec = roc_MASH_smooth$specificities,  stringsAsFactors = F)


# HCV Subclass
load("hcc_HCV_res.RData")
hcv_predictions <- lapply(TestPerformance.list, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "HCV",1,0))%>%mutate(Class2 = as.numeric(Class2)))
for(i in 1:length(hcv_predictions)) { hcv_predictions[[i]]$Index <- paste0("Model_",i) }
ROCs.computed <- lapply(hcv_predictions, function(x) with(x, roc(Class2 ~ One)))
ROCs.computed <- lapply(ROCs.computed, function(x) smooth(x, method = "density"))
ROCs.computed <- lapply(ROCs.computed, function(x) data.frame(Sens = x$sensitivities, Spec = x$specificities, stringsAsFactors = F))
AUCs <- lapply(hcv_predictions, function(x) with(x, roc(Class2 ~ One)$auc))
AUCs <- unlist(AUCs)
# Obtain AUROC
t.test(AUCs, conf.level = 0.95)
# Combine all model predictions into a large data frame
combined_predictions <- do.call(rbind, hcv_predictions)
# Compute the overall ROC
roc_hcv <- with(combined_predictions, roc(Class2 ~ One))
# Compute the smoothed ROC
roc_hcv_smooth <- smooth(roc_hcv, method = "density")
# Extract the sensitivities and specificities
roc_hcv <- data.frame(Sens = roc_hcv_smooth$sensitivities, Spec = roc_hcv_smooth$specificities, stringsAsFactors = F)

# ND Subclass
load("hcc_ND_res.RData")
nd_predictions <- lapply(TestPerformance.list, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "ND",1,0))%>%mutate(Class2 = as.numeric(Class2)))
for(i in 1:length(nd_predictions)) { nd_predictions[[i]]$Index <- paste0("Model_",i) }
ROCs.computed <- lapply(nd_predictions, function(x) with(x, roc(Class2 ~ One)))
ROCs.computed <- lapply(ROCs.computed, function(x) smooth(x, method = "density"))
ROCs.computed <- lapply(ROCs.computed, function(x) data.frame(Sens = x$sensitivities, Spec = x$specificities, stringsAsFactors = F))
AUCs <- lapply(nd_predictions, function(x) with(x, roc(Class2 ~ One)$auc))
AUCs <- unlist(AUCs)
# Obtain AUROC
t.test(AUCs, conf.level = 0.95)
# Combine all model predictions into a large data frame
combined_predictions <- do.call(rbind, nd_predictions)
# Compute the overall ROC
roc_nd <- with(combined_predictions, roc(Class2 ~ One))
# Compute the smoothed ROC
roc_nd_smooth <- smooth(roc_nd, method = "density")
# Extract the sensitivities and specificities
roc_nd <- data.frame(Sens = roc_nd_smooth$sensitivities, Spec = roc_nd_smooth$specificities,  stringsAsFactors = F)

# combine all ROC curves
roc_hbv$Dataset <- "HBV"
roc_ald$Dataset <- "ALD"
roc_MASH$Dataset <- "MASH"
roc_hcv$Dataset <- "HCV"
roc_nd$Dataset <- "ND"
roc_ctl$Dataset <- "CTL"

roc_all <- rbind(roc_hbv, roc_ald, roc_MASH, roc_hcv, roc_nd, roc_ctl)

# plot Fig.S3A ################################
ggplot(data = roc_all, aes(x = 1 - Spec, y = Sens, colour = Dataset)) +
 theme_bw() +
 geom_line(size = 5) +
 scale_color_manual(values = c("HBV" = "#d8222b", "ALD" = "#00aacb", "HCV" = "#fa9d48", "MASH" = "#424242", "ND" = "grey", "CTL" = "#5d6a85")) +
 theme(legend.position = "none",
       panel.border = element_rect(linewidth = 4, color = "black"),
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank(),
       axis.text.x = element_text(vjust = 0.5, color = "black"),  # Change the color here
       axis.text.y = element_text(hjust = -0.5, color = "black"),  # And here
       axis.ticks.length = unit(0.5, "cm")) +
 geom_abline(intercept = 0, slope = 1, color ="black", linetype = "dashed", size =1) +
 xlab("1 - Specificity") + 
 ylab("Sensitivity") +
 theme(text=element_text(size=70,color="black"),
       axis.title = element_text(size = 70, color="black"),
       plot.title = element_text(face = "italic", color="black", size=40),
       plot.margin = unit(c(0.3, 1, 0.3, 0.3), "cm"))

## PCA Analysis
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)

# load lcpm matrix
count <- readRDS("n236_lcpm.rds")

# Extract HBV DMR
load("hcc_HBV_res.RData")
## Select feature with threshold ############################
nFreq <- 1
Features_vec <- c()
for (j in 1:100) {
 kFold.list <- All.kFold[[j]]
 for (i in 1:10) {
  Features_vec <- c(Features_vec, kFold.list[[i]]$Model$finalModel$xNames)
 }
}
re <- as.data.frame(table(Features_vec))
df_temp <- re[re$Freq >= nFreq,]
HBV_DMR <- as.data.frame(df_temp$Features_vec)
dim(HBV_DMR)

# Extract HCV DMR
load("hcc_HCV_res.RData")
## Select feature with threshold ############################
nFreq <- 1
Features_vec <- c()
for (j in 1:100) {
 kFold.list <- All.kFold[[j]]
 for (i in 1:10) {
  Features_vec <- c(Features_vec, kFold.list[[i]]$Model$finalModel$xNames)
 }
}
re <- as.data.frame(table(Features_vec))
df_temp <- re[re$Freq >= nFreq,]
HCV_DMR <- as.data.frame(df_temp$Features_vec)
dim(HCV_DMR)

# Extract ALD DMR
load("hcc_ALD_res.RData")
## Select feature with threshold ############################
nFreq <- 1
Features_vec <- c()
for (j in 1:100) {
 kFold.list <- All.kFold[[j]]
 for (i in 1:10) {
  Features_vec <- c(Features_vec, kFold.list[[i]]$Model$finalModel$xNames)
 }
}
re <- as.data.frame(table(Features_vec))
df_temp <- re[re$Freq >= nFreq,]
ALD_DMR <- as.data.frame(df_temp$Features_vec)
dim(ALD_DMR)

# Extract MASH DMR
load("hcc_MASH_res.RData")
## Select feature with threshold ############################
nFreq <- 1
Features_vec <- c()
for (j in 1:100) {
 kFold.list <- All.kFold[[j]]
 for (i in 1:10) {
  Features_vec <- c(Features_vec, kFold.list[[i]]$Model$finalModel$xNames)
 }
}
re <- as.data.frame(table(Features_vec))
df_temp <- re[re$Freq >= nFreq,]
MASH_DMR <- as.data.frame(df_temp$Features_vec)
dim(MASH_DMR)

# Extract ND (non-liver disease) DMR
load("hcc_ND_res.RData")
## Select feature with threshold ############################
nFreq <- 1
Features_vec <- c()
for (j in 1:100) {
 kFold.list <- All.kFold[[j]]
 for (i in 1:10) {
  Features_vec <- c(Features_vec, kFold.list[[i]]$Model$finalModel$xNames)
 }
}
re <- as.data.frame(table(Features_vec))
df_temp <- re[re$Freq >= nFreq,]
# Final selected signature
ND_DMR <- as.data.frame(df_temp$Features_vec)
dim(ND_DMR)

# Extract CTL (cancer-free healthy controls) DMR
load("hcc_CTL_res.RData")
## Select feature with threshold ############################
nFreq <- 1
Features_vec <- c()
for (j in 1:100) {
 kFold.list <- All.kFold[[j]]
 for (i in 1:10) {
  Features_vec <- c(Features_vec, kFold.list[[i]]$Model$finalModel$xNames)
 }
}
re <- as.data.frame(table(Features_vec))
df_temp <- re[re$Freq >= nFreq,]
CTL_DMR <- as.data.frame(df_temp$Features_vec)
dim(CTL_DMR)

Features <- rbind(HBV_DMR,HCV_DMR,ALD_DMR,MASH_DMR,ND_DMR,CTL_DMR)

# Identify and remove all duplicated rows, including the first occurrence
all_duplicates <- duplicated(Features) | duplicated(Features, fromLast = TRUE)
# Print the number of duplicate rows to be removed
num_all_duplicates <- sum(all_duplicates)
print(paste("Number of rows to be removed (including the first occurrence):", num_all_duplicates))

# Remove all duplicated rows
Features_no_duplicates <- Features[!all_duplicates, ]
Features_no_duplicates <- as.data.frame(Features_no_duplicates)
# Print the number of rows before and after removing duplicates
print(paste("Number of rows before removing duplicates:", nrow(Features)))
print(paste("Number of rows after removing duplicates:", nrow(Features_no_duplicates)))

### PCA analysis
colnames(Features_no_duplicates) <- c("DMR")
DMR_count <- count[rownames(count) %in% Features_no_duplicates$DMR,]
sample <- readRDS('sample_124.rds')
pheno <- sample
table(pheno$Disease)
table(pheno$Cirr_group)
DMR_count <- DMR_count[,match(pheno$Sample_ID, colnames(DMR_count))]
identical(colnames(DMR_count), pheno$Sample_ID)

dmr <- as.data.frame(DMR_count)
dmr$bin <- rownames(dmr)
dmr <- dmr[,c(125,1:124)]
data <- as.data.frame(t(dmr[,c(2:125)]))
data$sample <- rownames(data)
data <- cbind(pheno, data[ ,c(25643,1:25642)])
identical(data$Sample_ID, data$sample)
res.pca <- PCA(data[ ,c(9:25650)], graph = FALSE)

# PCA plot
basic_plot <- fviz_pca_ind(res.pca, label="none")
b_d <- basic_plot$data
b_d <- cbind(b_d, data$Disease, data$Cirr_group,data$Sample_ID)
colnames(b_d)[c(7,8,9)] <- c("Group","Group_2","Sample")
b_d <- b_d[c(7,8,9,2,3)]

# show PC-1
print(get_eigenvalue(res.pca)[1,2])
# show PC-2
print(get_eigenvalue(res.pca)[2,2])

# PCA plot (Fig.S3B) ################################

qplot(data = b_d, x = x, y = y,  size = I(0), alpha = I(0.1)) +
 theme_bw() +
 theme(text = element_text(size = 65), axis.line = element_line(size = 1)) +
 theme(panel.border = element_rect(size = 2, color = "black"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       axis.ticks.length = unit(0.5, "cm")) +
 #panel.background = element_rect(fill = I("grey"))) +
 geom_point(data = b_d %>% filter(Group == "HBV" & Group_2 == "NC"), aes(y = y, x = x), shape = 21, size = 16, alpha = I(0.8), color = "black", fill = I("#d8222b"), stroke = 1) +
 geom_point(data = b_d %>% filter(Group == "HBV" & Group_2 == "Cirrhosis"), aes(y = y, x = x), shape = 21, size = 16, alpha = I(0.8), color = "black", fill = I("#d8222b"), stroke = 1) +
 geom_point(data = b_d %>% filter(Group == "HBV" & Group_2 == "Cirrhosis"), aes(y = y, x = x), shape = 21, size = 4, alpha = I(1), color = "black", fill = I("black"), stroke = 1.5) +
 geom_point(data = b_d %>% filter(Group == "ND" & Group_2 == "NC"), aes(y = y, x = x), shape = 21, size = 16, alpha = I(0.8), color = "black", fill = I("#c0c0c0"), stroke = 1) +
 geom_point(data = b_d %>% filter(Group == "ND" & Group_2 == "Cirrhosis"), aes(y = y, x = x), shape = 21, size = 16, alpha = I(0.8), color = "black", fill = I("#c0c0c0"), stroke = 1) +
 geom_point(data = b_d %>% filter(Group == "ND" & Group_2 == "Cirrhosis"), aes(y = y, x = x), shape = 21, size = 4, alpha = I(1), color = "black", fill = I("black"), stroke = 1.5) +
 geom_point(data = b_d %>% filter(Group == "MASH"), aes(y = y, x = x), shape = 21, size = 15, alpha = I(0.9), color = "black", fill = I("#545454"), stroke = 1) +
 geom_point(data = b_d %>% filter(Group == "ALD" & Group_2 == "NC"), aes(y = y, x = x), shape = 21, size = 16, alpha = I(0.9), color = "black", fill = I("#00aacb"), stroke = 1) +
 geom_point(data = b_d %>% filter(Group == "ALD" & Group_2 == "Cirrhosis"), aes(y = y, x = x), shape = 21, size = 16, alpha = I(0.9), color = "black", fill = I("#00aacb"), stroke = 1) +
 geom_point(data = b_d %>% filter(Group == "ALD" & Group_2 == "Cirrhosis"), aes(y = y, x = x), shape = 21, size = 4, alpha = I(1), color = "black", fill = I("black"), stroke = 1.5) +
 geom_point(data = b_d %>% filter(Group == "HCV" & Group_2 == "NC"), aes(y = y, x = x), shape = 21, size = 16, alpha = I(0.9), color = "black", fill = I("#fb9e49"), stroke = 1) +
 geom_point(data = b_d %>% filter(Group == "HCV" & Group_2 == "Cirrhosis"), aes(y = y, x = x), shape = 21, size = 16, alpha = I(0.9), color = "black", fill = I("#fb9e49"), stroke = 1) +
 geom_point(data = b_d %>% filter(Group == "HCV" & Group_2 == "Cirrhosis"), aes(y = y, x = x), shape = 21, size = 4, alpha = I(1), color = "black", fill = I("black"), stroke = 1.5) +
 geom_point(data = b_d %>% filter(Group == "CTL"), aes(y = y, x = x), shape = 0, size = 12, alpha = I(0.8), color = "#0a234e", fill = I("#0a234e"), stroke = 2) +
 geom_point(data = b_d %>% filter(Group == "Autoimmune" & Group_2 == "NC"), aes(y = y, x = x), shape = 5, size = 12, alpha = I(1), color = "blue", fill = I("#942193"), stroke = 2) +
 geom_point(data = b_d %>% filter(Group == "Autoimmune" & Group_2 == "Cirrhosis"), aes(y = y, x = x), shape = 9, size = 12, alpha = I(1), color = "blue", fill = I("#942193"), stroke = 2) +
 geom_point(data = b_d %>% filter(Group == "PBC"), aes(y = y, x = x), shape = 14, size = 12, alpha = I(1), color = "#942193", fill = I("#0051ca"), stroke = 2) +
 geom_point(data = b_d %>% filter(Group == "HBV-MASH" & Group_2 == "NC"), aes(y = y, x = x), shape = 24, size = 12, alpha = I(0.6), color = "black", fill = I("#963b40"), stroke = 1) +
 geom_point(data = b_d %>% filter(Group == "HBV-MASH" & Group_2 == "Cirrhosis"), aes(y = y, x = x), shape = 24, size = 12, alpha = I(0.6), color = "black", fill = I("#963b40"), stroke = 1) +
 geom_point(data = b_d %>% filter(Group == "HBV-MASH" & Group_2 == "Cirrhosis"), aes(y = y, x = x), shape = 21, size = 4, alpha = I(1), color = "black", fill = I("black"), stroke = 1.5) +
 geom_point(data = b_d %>% filter(Group == "HBV-ALD"), aes(y = y, x = x), shape = 22, size = 14, alpha = I(0.6), color = "black", fill = I("#942193"), stroke = 2) +
 geom_point(data = b_d %>% filter(Group == "HBV-ALD"), aes(y = y, x = x), shape = 21, size = 4, alpha = I(1), color = "black", fill = I("black"), stroke = 1.5) +
 geom_point(data = b_d %>% filter(Group == "HCV-ALD"), aes(y = y, x = x), shape = 23, size = 14, alpha = I(0.7), color = "black", fill = I("grey"), stroke = 1) +
 geom_point(data = b_d %>% filter(Group == "HCV-ALD"), aes(y = y, x = x), shape = 21, size = 4, alpha = I(1), color = "black", fill = I("black"), stroke = 1.5) +
 geom_point(data = b_d %>% filter(Group == "MASH-ALD"), aes(y = y, x = x), shape = 25, size = 12, alpha = I(0.4), color = "black", fill = I("#10a1a1"), stroke = 1) +
 geom_point(data = b_d %>% filter(Group == "MASH-ALD"), aes(y = y, x = x), shape = 6, size = 12, alpha = I(0.4), color = "black", fill = I("black"), stroke = 1) +
 xlab("PC1 (7.80%)")+ ylab("PC2 (5.26%)")+
 theme(plot.margin = unit(c(0.3, 2.5, 0.5, 0.3), "cm")) +
 theme(axis.text.x = element_text(size = 55, color = "black"), axis.text.y = element_text(size = 55, color = "black"),
       axis.title.x = element_text(vjust = 0, size = 55, color = "black"), axis.title.y = element_text(vjust = -1.5, size = 55, color = "black"))
