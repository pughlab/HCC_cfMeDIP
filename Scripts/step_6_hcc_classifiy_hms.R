# Description: This script bulid the HCC classifier using kFold cross-validation machine learning strategy.
## Author：Kui Chen → kui.chen@uhn.ca; Ping Luo → ping.luo@uhn.ca #############################################
## Date Created: 2022-07-18 ###################################################################################
## Last Modified: 2024-02-10 ##################################################################################
# Runnning ~10 hours with RStudio at macOS with 64GB memory configuration #####################################
# Optional : Runnning ~2 days at H4H cluster with ≥ 180GB memory (veryhimem) configuration using R/3.5.0 ######
# Before starting, do the followings ##########################################################################
# 1.make working directory "/path/to/your/ML" #################################################################
# 2.creat soft link of "n236_cpm.rds", "n236_lcpm.rds" in the working directory ###############################
# 3.copy support files "sample_236.txt", "black_bin_v2.RData" to the working directory ########################
###############################################################################################################

rm(list=ls())
gc()
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(caret)
library(randomForest)
library(glmnet)
library(limma)
library(doParallel)
Cluster <-makeCluster(24)
Cluster <- registerDoParallel(Cluster)

####################################################################
#### kFold cross-validation machine learning #######################
####################################################################

setwd("/path/to/your/counts")

# Prepare matrix
load("n236_lcpm.RData")

# exclude validation cohort data
sample <- sample[which(!sample$Group_2 == "V"),]
table(sample$Group_1)
table(sample$Group_2)
table(sample$Group_3)
table(sample$Group_4)

matrix <- matrix[ ,match(sample$RUN_ID,colnames(matrix))]
identical(colnames(matrix), sample$RUN_ID)
dim(matrix)

#filter out 300-bp bins located in the ENCODE Blacklist (Problematic Regions of the Human Genome)
load("black_bin_v2.RData")
matrix <- subset(matrix, !rownames(matrix) %in% black_bin$black_bin)
identical(colnames(matrix), sample$RUN_ID)
dim(matrix)

# setup groups
hcc <- matrix[, sample$Group_1 == 'HCC']
ctr <- matrix[, sample$Group_1 == 'CTL']
rm(count)
dim(hcc)
dim(ctr)

# 10 fold cross-validation modeling

Feature_list <- list()
Prob_list <- list()
Y_true_list <- list()
Group_list <- list()

Y_sub <-  c(rep("hcc", ncol(hcc)), rep("ctr", ncol(ctr)))

numF <- 150
numK <- 100

for (k in 1:numK) {

 Data_sub <- cbind(hcc, ctr)
 Prob <- NULL
 Y_true <- NULL
 Features <- NULL
 Groups <- NULL
 Y_pred <- NULL

 samples <- createFolds(Y_sub, k = 10, list = T, returnTrain=T)

 # 10 fold CV
 for (j in 1:10) {
  X_train <- Data_sub[,samples[[j]]]
  y_train <- Y_sub[samples[[j]]]

  test_ind <- setdiff(seq(1, length(Y_sub)), samples[[j]])

  X_test <- Data_sub[,test_ind]
  y_test <- Y_sub[test_ind]

  # Find DMR
  DMR.classes <- y_train
  Des <- model.matrix(~0 + DMR.classes)
  colnames(Des) <- levels(factor(DMR.classes))

  LimmaFit <- lmFit(X_train, Des)%>%
   contrasts.fit(., makeContrasts(hcc-ctr, levels = Des))%>%
   eBayes(., trend = TRUE)%>%
   topTable(., number = nrow(X_train))
  LimmaFit <- LimmaFit%>%.[order(.$t,decreasing = T),]

  TotalRows <- nrow(LimmaFit) - (numF - 1)
  Feature_vec <- rownames(rbind(LimmaFit[1:numF,],LimmaFit[TotalRows:nrow(LimmaFit),]))
  Group_vec <- c(rep('hcc', numF), rep('ctr', numF))

  # Model training and prediction
  ind <- match(Feature_vec, rownames(Data_sub))
  mtx <- X_train[ind,]
  mtx <- as.data.frame(t(mtx))
  mtx$Type <- y_train
  mtx[is.na(mtx)] <- 0

  control <- trainControl(method='repeatedcv', number=10, repeats=1)
  tunegrid <- expand.grid(.mtry=c(10))
  Model <- train(Type~., data=mtx,
                 method='rf',
                 metric='Accuracy',
                 tuneGrid=tunegrid,
                 trControl=control)

  prob <- predict(Model, newdata = t(X_test[ind,]), type = "prob") %>% data.frame
  y_pred <- predict(Model, newdata = t(X_test[ind,]), type = "raw") %>% as.vector

  message(paste0(length(X_train)," x Train"))
  message(colnames(X_train))
  message(paste0(length(X_test)," x Test"))
  message(colnames(X_test))
  message(paste0(length(Feature_vec), "s DMR"))
  message(paste0("cycle ", k , " NO.", j ," of DMR selection done"))

  # save result
  if (j == 1) {
   Prob <- prob
   Y_true <- y_test
   Features <- Feature_vec
   Groups <- Group_vec
   Y_pred <- y_pred
  } else {
   Prob <- rbind(Prob, prob)
   Y_true <- c(Y_true, y_test)
   Features <- c(Features, Feature_vec)
   Groups <- c(Groups, Group_vec)
   Y_pred <- c(Y_pred, y_pred)
  }
 }

 Prob_list[[k]] <- Prob
 Y_true_list[[k]] <- Y_true
 Feature_list[[k]] <- Features
 Group_list[[k]] <- Groups

 Prob_list[[k]]$ActualClass <- Y_true
 Prob_list[[k]]$PredictedClass <- Y_pred
}

save(Prob_list, Y_true_list, Feature_list, Group_list, sample, file = "kFold_ML_res.RData")

##########################################################################################
#### Evaluate modeling performance with a validation cohort ############################
##########################################################################################

# get Validation cohort HMS using slected DMRs list
load("kFold_ML_res.RData")

# organize matrix
count <- readRDS("n236_lcpm.rds")
count <- as.data.frame(count)
sample <- read.table(file="sample_236k.txt", sep = "\t", header=T, stringsAsFactors=FALSE)
dim(count)

#check how many DMR in black_list
load("black_bin_v2.RData")
count <- subset(count, !rownames(count) %in% black_bin$black_bin)
dim(count)

# training-testing groups setup
matrix <- count[ ,match(sample$Sample_ID,colnames(count))]
identical(colnames(matrix), sample$Sample_ID)
table(sample$Group_1)
table(sample$Group_2)
table(sample$Group_3)
table(sample$Group_4)

# Setup training and testing dataset groups
Matrix.train.all <- matrix[,sample$Group_2 == "D"]
Matrix.test.all <- matrix[,sample$Group_2 == "V"]

# Generating train phenotype
Phenotype.train <- sample[sample$Group_2 == "D",]
df.train <- data.frame(ID = colnames(Matrix.train.all), Classes = Phenotype.train$Group_1)
Features.CVparam<- trainControl(method = "repeatedcv",number = 10, repeats =3, verboseIter=TRUE,returnData=FALSE, classProbs = TRUE,savePredictions=FALSE)
NewAnn <- ifelse(df.train$Classes == "HCC", "One","Others")
Features <- unique(unlist(Feature_list))
Features <- unlist(Feature_list)
mtry.val <- length(Features)/10

Val.prob.list <- list()

for (i in 1:100) {
 Features <- Feature_list[[i]]

 # Generating train and test DMRs matrix
 Matrix.train <- Matrix.train.all[Features,]
 Matrix.test <- Matrix.test.all[Features,]

 # Training the model
 Model <- train(x = t(Matrix.train), y = factor(NewAnn), trControl = Features.CVparam, method = "rf" , tuneGrid = expand.grid(.mtry = mtry.val), metric = "Kappa")
 message(i)
 # Start prediction
 Prediction.classProbs <- predict(Model, newdata = t(Matrix.test), type = "prob") %>% data.frame
 Prediction.classProbs$ActualClass <- sample[sample$Group_2 %in% c("V"),]$Group_1
 Prediction.classProbs$SubClass <- sample[sample$Group_2 %in% c("V"),]$Group_3
 Prediction.classProbs$PredictedClass <- predict(Model, newdata = t(Matrix.test), type = "raw")
 Val.prob.list[[i]] <- list(Model = Model, TestPred = Prediction.classProbs)
}

saveRDS(Val.prob.list, 'Val.prob.list.rds')

Val.prob.list <- readRDS('Val.prob.list.rds')
# generate HMS of Validation cohort
Val.prob.list <- lapply(Val.prob.list, "[[", "TestPred")
Val_score_list <- list()
for (i in 1:100)
{
 Val_score_list[[i]] <- Val.prob.list[[i]]
 Val_score_list[[i]]$Sample_ID <- rownames(Val_score_list[[i]])
 Val_score_list[[i]] <- Val_score_list[[i]][order(Val_score_list[[i]]$Sample_ID),]
 Val_score_list[[i]] <- Val_score_list[[i]]$One
}

# attach #1 Prob
Prob_1 <- as.data.frame(Val.prob.list[[1]])
Prob_1$Sample_ID <- rownames(Prob_1)
Prob_1 <- Prob_1[order(Prob_1$Sample_ID),]
V_score <- cbind(Prob_1,as.data.frame(Val_score_list))
V_score[1:10,1:10]
# double-check the cbind by confirm 1st one
identical(V_score[,1],V_score[,7])
identical(rownames(V_score),V_score$Sample_ID)

# generate score data.frame
V_score <- V_score[,-c(1:6)]
colnames(V_score) <- paste0("Prob_",1:100)
dim(V_score)

# get mean score
V_score$Mean <- apply(V_score[,1:100],1,mean)
V_score$Sample <- rownames(V_score)
V_score <- V_score[,c(102,101,1:100)]


# # get Discovery cohort HMS
D_score_list <- list()
for (i in 1:100) {
 D_score_list[[i]] <- Prob_list[[i]]
 D_score_list[[i]]$Sample_ID <- rownames(D_score_list[[i]])
 D_score_list[[i]] <- D_score_list[[i]][order(D_score_list[[i]]$Sample_ID),]
 D_score_list[[i]] <- D_score_list[[i]]$hcc
}

# attach #1 Prob
Prob_1 <- as.data.frame(Prob_list[[1]])
Prob_1$Sample_ID <- rownames(Prob_1)
Prob_1 <- Prob_1[order(Prob_1$Sample_ID),]
D_score <- cbind(Prob_1,as.data.frame(D_score_list))

# double-check the cbind by confirm 1st one
identical(D_score[,2],D_score[,6])

# generate score data.frame
D_score <- D_score[,-c(1:5)]
colnames(D_score) <- paste0("Prob_",1:100)
dim(D_score)

# get mean score
D_score$Mean <- apply(D_score[,1:100],1,mean)
D_score$Sample <- rownames(D_score)
D_score <- D_score[,c(102,101,1:100)]

# Save HCC Methylation Score (HMS) Results
## Fig.2F, Fig.2G, Fig.3, Fig.4C, Fig.5 and Fig.S3C were generated based on the HMS results ########################
sample <- read.table(file="sample_236k.txt", sep = "\t", header=T, stringsAsFactors=FALSE)
sq <- read.table(file="sample_sq.txt", sep = "\t", header=T, stringsAsFactors=FALSE)

All_score <- rbind(D_score, V_score)
HMS <- merge(sample, All_score, by.x = "Sample_ID", by.y = "Sample", all = TRUE)
names(HMS)[names(HMS) == "Mean"] <- "HMS"
D_score <- HMS[HMS$Group_2 == "D", ]
V_score <- HMS[HMS$Group_2 == "V", ]
HMS_sq <- HMS[HMS$Sample_ID %in% sq$Sample_Sq, ]

library(openxlsx)
# Create a new Workbook object
wb <- createWorkbook()

# Add a sheet for each data frame
addWorksheet(wb, "HMS")
addWorksheet(wb, "HMS_sq")
addWorksheet(wb, "D_score")
addWorksheet(wb, "V_score")
# Write each data frame to its corresponding sheet
writeData(wb, sheet = "HMS", HMS)
writeData(wb, sheet = "HMS_sq", HMS_sq)
writeData(wb, sheet = "D_score", D_score)
writeData(wb, sheet = "V_score", V_score)
# Save the Workbook to disk
saveWorkbook(wb, "Table_S5_HMS.xlsx", overwrite = TRUE)

# Obtain sensitivity, specificity, and accuracy of modeling performance of discovery cohort samples
## data for generating Fig.2E (top) #############################
library(caret)
load("kFold_ML_res.RData")
HCCs.predictions <- lapply(Prob_list, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "hcc",1,0))%>%mutate(Class2 = as.numeric(Class2)))
for(i in 1:length(HCCs.predictions)) {  HCCs.predictions[[i]]$Index <- paste0("Model_",i)}

confusion.matrix.list <- lapply(1:length(HCCs.predictions), function(i) {
 confusionMatrix(reference = factor(HCCs.predictions[[i]]$Class2), data = factor(ifelse(HCCs.predictions[[i]]$PredictedClass == "hcc",1,0)))
})

confusion.matrix <- lapply(confusion.matrix.list, function(x) x$table)
confusion.matrix <- Reduce(function(x, y) x + y, result)
confusion.matrix

# Obtain the accuracy of f-HCC samples predicted as HCC ###############
# data for generating Fig.2E (bottom) #############################
Val.prob.list <- readRDS('Val.prob.list.rds')

# Obtain accuracy of all follow-up HCC (f-HCC) samples predicted as HCC
hcc_samples <- sum(sapply(Val.prob.list, function(x) sum(x$TestPred$SubClass == "V-HCC")))
hcc_correct_predictions <- sum(sapply(Val.prob.list, function(x) sum(x$TestPred$SubClass == "V-HCC" & x$TestPred$PredictedClass == "One")))
print(paste("f-HCCs predicted as HCC:", hcc_correct_predictions))
print(paste("Total predictions:", hcc_samples))
if(hcc_samples > 0) {
 hcc_accuracy <- hcc_correct_predictions / hcc_samples
 print(paste("Accuracy of all f-HCC samples were predicted as HCC:", hcc_accuracy))
} else {
 print("No HCC samples detected.")
}

# Obtain accuracy of recurrence f-HCC samples predicted as HCC
rec_samples <- sum(sapply(Val.prob.list, function(x) sum(x$TestPred$SubClass == "Rk")))
rec_correct_predictions <- sum(sapply(Val.prob.list, function(x) sum(x$TestPred$SubClass == "Rk" & x$TestPred$PredictedClass == "One")))
print(paste("Recurrence f-HCCs predicted as HCC :", rec_correct_predictions))
print(paste("Total predictions:", rec_samples))
if(hcc_samples > 0) {
 rec_accuracy <- rec_correct_predictions / rec_samples
 print(paste("Accuracy of Recurrence f-HCC samples were predicted as HCC:", rec_accuracy))
} else {
 print("No HCC samples detected.")
}

# Obtain accuracy of remission f-HCC samples predicted as HCC
rem_samples <- sum(sapply(Val.prob.list, function(x) sum(x$TestPred$SubClass == "Rm")))
rem_correct_predictions <- sum(sapply(Val.prob.list, function(x) sum(x$TestPred$SubClass == "Rm" & x$TestPred$PredictedClass == "One")))
print(paste("Remission f-HCCs predicted as HCC:", rem_correct_predictions))
print(paste("Total predictions:", rem_samples))
if(hcc_samples > 0) {
 rem_accuracy <- rem_correct_predictions / rem_samples
 print(paste("Accuracy of Remission f-HCC samples were predicted as HCC:", rem_accuracy))
} else {
 print("No HCC samples detected.")
}


################################################
#### AUROC #####################################
################################################
# obtain mean AUROC
library(pROC)
#load("kFold_ML_res.RData")
# method 1
auc_list <- lapply(1:length(Prob_list), function(i) {
 roc_obj <- roc(Y_true_list[[i]], Prob_list[[i]]$hcc)
 auc(roc_obj)
})

mean_auc <- mean(unlist(auc_list))
auc <- unlist(auc_list)

alpha <- 0.05
t_critical <- qt(1 - alpha / 2, df = length(auc) - 1)
se_auc <- sd(auc) / sqrt(length(auc))
lower_bound <- mean_auc - t_critical * se_auc
upper_bound <- mean_auc + t_critical * se_auc
Mean_AUC <- paste("Mean AUC:", round(mean_auc, 3), "95% CI:", round(lower_bound, 3), "-", round(upper_bound, 3))
Mean_AUC

# method 2
## Get 100 individual curves data (Mean AUROC)
HCCs.predictions <- lapply(Prob_list, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "hcc",1,0))%>%mutate(Class2 = as.numeric(Class2)))
for(i in 1:length(HCCs.predictions)) {  HCCs.predictions[[i]]$Index <- paste0("Model_",i)                    }
ROCs.computed <- lapply(HCCs.predictions, function(x) with(x, roc(Class2 ~ hcc)))
ROCs.computed <- lapply(ROCs.computed, function(x) smooth(x, method = "density"))
ROCs.computed <- lapply(ROCs.computed, function(x) data.frame(Sens = x$sensitivities, Spec = x$specificities, stringsAsFactors = F))
AUCs <- lapply(HCCs.predictions, function(x) with(x, roc(Class2 ~ hcc)$auc))
AUCs <- unlist(AUCs)
t.test(AUCs, conf.level = 0.95)

## generate ROC curve (Fig.2A) ##########
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(grid)
pal <- colorRampPalette(c("black","#db046e"))(100)
pal <- colorRampPalette(c("#526c9c","#671443"))(100)
pal <- colorRampPalette(c("black","red"))(100)
for (i in 1:length(ROCs.computed)) { ROCs.computed[[i]]$Index <- paste0("Model_",i) ; ROCs.computed[[i]]$auc <- AUCs[[i]] }
ROCs.computed.bound <- do.call(rbind, ROCs.computed)
ROCs.computed.bound <- arrange(ROCs.computed.bound, auc)%>%
 mutate(Index = factor(Index, ordered = TRUE , levels = unique(Index)))

ggplot(data = ROCs.computed.bound, aes(x = 1 - Spec, y = Sens, colour = Index , label = NULL))+
 theme_bw() +
 geom_step(alpha = I(1), size = I(1))+
 theme(legend.position = "none",
       panel.border = element_rect(linewidth = 4, color = "black"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       axis.text.x = element_text(vjust = 1, color = "black"),  # Increase gap on x-axis
       axis.text.y = element_text(hjust = -0.5,color = "black"),  # Increase gap on y-axis
       axis.ticks.length = unit(0.5, "cm"))+
 geom_abline(intercept = 0, slope = 1, color ="grey", linetype = "solid", size =2.5) +
 scale_colour_manual(values = pal)+
 xlab("1 - Specificity")+ 
 ylab("Sensitivity")+
 theme(text=element_text(size=70,color="black"),
       axis.title = element_text(size = 70, color="black"),
       plot.title = element_text(color="black", size=40),
       plot.margin = unit(c(0.3, 1, 0.3, 0.3), "cm")) +
       #opbtian the exact mean AUROC & 95%CI data from the "obtain mean AUROC" section
 annotate("text", label = "Mean AUROC = 0.999", x = 0.63, y = 0.14, size = 22, colour = "black") +
 annotate("text", label = "(95%CI 0.998 ~ 1)", x = 0.7, y = 0.04, size = 22,  colour = "black")


##########################################
#### PCA Analysis ########################
##########################################

library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)

load("kFold_ML_res.RData")

# Isolated DMRs
## get DMR genome coordinators
DMR_Coo_list <- lapply(1:100, function(x) {
 data.frame(Group = Group_list[[x]], Feature = Feature_list[[x]])
})

DMR_Coo <- do.call(rbind, DMR_Coo_list)

## confirm that all DMR are just assigened in single group (either hcc or ctr)
## Group the DMR_Coo dataframe by Feature and count the number of unique groups for each feature
group_count <- DMR_Coo %>% group_by(Feature) %>% summarize(n_groups = n_distinct(Group))

## Check if all features belong to only one group
if (all(group_count$n_groups == 1)) {
 print("All features belong to only one group")
} else if (any(group_count$n_groups == 2)) {
 print("Some features belong to both groups")
} else {
 print("Error: unexpected n_groups values")
}

## how many DMRs identified
unique_features <- unique(DMR_Coo$Feature)
# Create a new data.frame with unique feature values
unique_DMR_Coo <- data.frame(Group = DMR_Coo$Group[match(unique_features, DMR_Coo$Feature)], Feature = unique_features)
# check the quantities of hyper / hypo methylated DMRs
table(unique_DMR_Coo$Group)
hyper_DMR <- unique_DMR_Coo[which(unique_DMR_Coo$Group == 'hcc'),]
hypo_DMR <- unique_DMR_Coo[which(unique_DMR_Coo$Group == 'ctr'),]
nrow(hyper_DMR)
nrow(hypo_DMR)

## isolate DMRs appears at 1 time during 100 times of training-testing iterations
nFreq <- 1
DMR_vec <- c()
for (i in 1:100) {
 DMR_vec <- c(DMR_vec, Feature_list[[i]])
}
re <- as.data.frame(table(DMR_vec))
range(re$Freq)
df_temp <- re[re$Freq >= nFreq,]
DMR_selected <- as.data.frame(df_temp$DMR_vec)
dim(DMR_selected)
colnames(DMR_selected) <- "DMR_selected"

## assign DMR group info
merged_df <- merge(unique_DMR_Coo, DMR_selected, by.x = "Feature", by.y = "DMR_selected")
colnames(merged_df) <- c('DMR','Group')
# extract hyper/hypo DMR
hyper <- merged_df[which(merged_df$Group == 'hcc'),]
hypo <- merged_df[which(merged_df$Group == 'ctr'),]

## check the quantities of total, hyper- or hypo-methylated DMRs
nrow(merged_df)
nrow(hyper)
nrow(hypo)

## load matrix & sample group data for PCA analysis
load("n236_lcpm.RData")
DMR_count <- matrix[rownames(matrix) %in% DMR_selected$DMR_selected,]
table(sample$Group_1)
table(sample$Group_2)
table(sample$Group_3)
table(sample$Group_4)

#### PCA Analysis on discovery cohort (52 b-HCC + 35 CTL) only (to generate Fig.2B) ############################################################################################################
pheno <- sample[which(sample$Group_2 == 'D'),]
dmr <- DMR_count[,match(pheno$RUN_ID, colnames(DMR_count))]
identical(colnames(dmr), pheno$RUN_ID)
dmr$bin <- rownames(dmr)
dmr <- dmr[,c(88,1:87)]
#colnames(dmr)
data <- as.data.frame(t(dmr[,c(2:88)]))
data$sample <- rownames(data)
data <- cbind(pheno, data[ ,c(4995,1:4994)])
identical(data$RUN_ID, data$sample)
res.pca <- PCA(data[ ,c(8:5001)], graph = F)

# generate scree plot with top100 individuals and contributions of variables
fviz_contrib(res.pca, choice="var", axes =1:2, top =100)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

# plot PCA clusters on discovery cohort
basic_plot <- fviz_pca_ind(res.pca, label="none")
b_d <- basic_plot$data
b_d <- cbind(b_d, data$Group_1, data$Group_2, data$Group_3,data$RUN_ID)
colnames(b_d)[c(7,8,9,10)] <- c("Group_1","Group_2", "Group_3","Sample")
b_d <- b_d[c(7,8,9,10,2,3)]

table(sample$Group_1)

# show PC-1
print(get_eigenvalue(res.pca)[1,2])
# show PC-2
print(get_eigenvalue(res.pca)[2,2])

# add PC-1/PC-2 value manually then plot Figure 2B
qplot(data=b_d, x=x, y=y, color = "#ffffff", size = I(12), alpha = I(0)) +
 theme_bw() +
 theme(text=element_text(size=65), axis.line = element_line(size = 1))+
 theme(panel.border = element_rect(size = 2, color = "black"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       axis.ticks.length = unit(0.5, "cm"))+
 geom_point(data = b_d%>%filter(Group_1 == "CTL"), aes(y = y, x = x), shape =21, size =12, alpha =I(0.8), color = "black", fill = I("#0a3161"), stroke = 2)+
 geom_point(data = b_d%>%filter(Group_1 == "HCC"), aes(y = y, x = x), shape =21, size = 12, alpha =I(1), color="black", fill = I("#f3b32a"), stroke = 2)+
 xlab("PC1 (25.71%)")+ ylab("PC2 (3.67%)")+
 theme(plot.margin = unit(c(0.5, 1, 0.5, 0.3), "cm")) +
 scale_color_manual(name = "Group", values = c("HCC" = "#f3b229", "CTL" = "#0a3161"))+
 guides(color = guide_legend(override.aes = list(size = 10,alpha = I(1))))+
 theme(axis.text.x = element_text(size = 55,color="black"),axis.text.y = element_text(size = 55,color="black"),
       axis.title.x = element_text(vjust=0, size = 55, color="black"), axis.title.y = element_text(vjust=-1, size = 55, color="black"))

#### PCA Analysis on 35 CTL and 37 b-HCC in validaiton cohort (generate Fig.2C) #############################################################################################################
pheno <- sample[which(sample$Group_3 == 'CTL' | sample$Group_3 == 'V-HCC'),]

dmr <- DMR_count[,match(pheno$RUN_ID, colnames(DMR_count))]
identical(colnames(dmr), pheno$RUN_ID)
dmr$bin <- rownames(dmr)
dmr <- dmr[,c(73,1:72)]

#colnames(dmr)
data <- as.data.frame(t(dmr[,c(2:73)]))
data$sample <- rownames(data)

data <- cbind(pheno, data[ ,c(4995,1:4994)])
identical(data$RUN_ID, data$sample)

res.pca <- PCA(data[ ,c(8:5001)], graph = F)

# generate scree plot with top100 individuals and contributions of variables
fviz_contrib(res.pca, choice="var", axes =1:2, top =100)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

# plot PCA clusters on discovery cohort
basic_plot <- fviz_pca_ind(res.pca, label="none")

b_d <- basic_plot$data
b_d <- cbind(b_d, data$Group_1, data$Group_2, data$Group_3, data$RUN_ID)
colnames(b_d)[c(7,8,9,10)] <- c("Group_1","Group_2", "Group_3","Sample")
b_d <- b_d[c(7,8,9,10,2,3)]

table(sample$Group_1)

# show PC-1
print(get_eigenvalue(res.pca)[1,2])
# show PC-2
print(get_eigenvalue(res.pca)[2,2])

# add PC-1/PC-2 value manually then plot Fig.2C
qplot(data=b_d, x=x, y=y, color = "#ffffff", size = I(12), alpha = I(0)) +
 theme_bw() +
 theme(text=element_text(size=65), axis.line = element_line(size = 1))+
 theme(panel.border = element_rect(size = 2, color = "black"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       axis.ticks.length = unit(0.5, "cm"))+
 geom_point(data = b_d%>%filter(Group_1 == "CTL"), aes(y = y, x = x), shape =21, size =12, alpha =I(0.8), color = "black", fill = I("#0a3161"), stroke = 2)+
 geom_point(data = b_d%>%filter(Group_1 == "HCC"), aes(y = y, x = x), shape =21, size = 12, alpha =I(1), color="black", fill = I("#f3b32a"), stroke = 2)+
 xlab("PC1 (20.87%)")+ ylab("PC2 (6.74%)")+
 theme(plot.margin = unit(c(0.5, 1, 0.5, 0.3), "cm")) +
 scale_color_manual(name = "Group", values = c("HCC" = "#f3b229", "CTL" = "#0a3161"))+
 guides(color = guide_legend(override.aes = list(size = 10,alpha = I(1))))+
 theme(axis.text.x = element_text(size = 55,color="black"),axis.text.y = element_text(size = 55,color="black"),
       axis.title.x = element_text(vjust=0, size = 55, color="black"), axis.title.y = element_text(vjust=-1, size = 55, color="black"))

#### PCA Analysis on all samples (89 b-HCC,35 CTL & 112 f-HCC (generate Fig.2D) #########################################################################################################
pheno <- sample
dmr <- DMR_count[,match(pheno$RUN_ID, colnames(DMR_count))]
identical(colnames(dmr), pheno$RUN_ID)
dmr$bin <- rownames(dmr)
dmr <- dmr[,c(237,1:236)]
#colnames(dmr)
data <- as.data.frame(t(dmr[,c(2:237)]))
data$sample <- rownames(data)
data <- cbind(pheno, data[ ,c(4995,1:4994)])
identical(data$RUN_ID, data$sample)
res.pca <- PCA(data[ ,c(8:5001)], graph = F)

# generate scree plot with top100 individuals and contributions of variables
fviz_contrib(res.pca, choice="var", axes =1:2, top =100)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

# plot PCA clusters on discovery cohort
basic_plot <- fviz_pca_ind(res.pca, label="none")
b_d <- basic_plot$data
b_d <- cbind(b_d, data$Group_1, data$Group_2, data$Group_3, data$RUN_ID)
colnames(b_d)[c(7,8,9,10)] <- c("Group_1","Group_2", "Group_3","Sample")
b_d <- b_d[c(7,8,9,10,2,3)]

table(b_d$Group_3)

# show PC-1
print(get_eigenvalue(res.pca)[1,2])
# show PC-2
print(get_eigenvalue(res.pca)[2,2])

# add PC-1/PC-2 value manually then plot Figure 2D
qplot(data = b_d, x = x, y = y, color = "#ffffff", size = I(9), alpha = I(0)) +
 theme_bw() +
 theme(text = element_text(size = 65), axis.line = element_line(size = 1)) +
 theme(panel.border = element_rect(size = 2, color = "black"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       axis.ticks.length = unit(0.5, "cm")) +
 geom_point(data = b_d %>% filter(Group_3 == "Rm"), aes(y = y, x = x), shape = 21, size = 10, alpha = I(1), color = "black", fill = I("#1f78b3"), stroke = 2) +
 geom_point(data = b_d %>% filter(Group_3 == "CTL"), aes(y = y, x = x), shape = 21, size = 10, alpha = I(0.2), color = "black", fill = I("#0a3161"), stroke = 2) +
 geom_point(data = b_d %>% filter(Group_3 == "CTL"), aes(y = y, x = x), shape = 1, size = 10, alpha = I(1), color = "black", fill = I("#0a3161"), stroke = 2) +
 geom_point(data = b_d %>% filter(Group_3 == "D-HCC"), aes(y = y, x = x), shape = 21, size = 10, alpha = I(0.9), color = "black", fill = I("#f3b229"), stroke = 2) +
 geom_point(data = b_d %>% filter(Group_3 == "V-HCC"), aes(y = y, x = x), shape = 21, size = 10, alpha = I(0.9), color = "black", fill = "#f3b229", stroke = 2) +
 geom_point(data = b_d %>% filter(Group_3 == "Rk"), aes(y = y, x = x), shape = 21, size = 10, alpha = I(0.9), color = "black", fill = I("#b31842"), stroke = 2) +
 xlab("PC1 (13.81%)")+ ylab("PC2 (5.32%)")+
 theme(plot.margin = unit(c(1, 1, 0.5, 0.3), "cm")) +
 scale_color_manual(name = "Group", values = c("CTL" = "#0a3161", "D-HCC" = "#f3b229", "V-HCC" = "#f3b229", "Rk" = "#b31842", "Rm" = "#1f78b3")) +
 theme(axis.text.x = element_text(size = 55, color = "black"), axis.text.y = element_text(size = 55, color = "black"),
       axis.title.x = element_text(vjust = 0, size = 55, color = "black"), axis.title.y = element_text(vjust = -0.5, size = 55, color = "black"))

