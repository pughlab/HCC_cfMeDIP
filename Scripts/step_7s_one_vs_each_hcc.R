# Description: This code supports the one-vs-each HCC subtyping main script
## Author：Ping Luo → ping.luo@uhn.ca, Kui Chen → kui.chen@uhn.ca
## Date Created: 2022-09-2
## Last Modified: 2023-12-16

SplitkFold <- function(Mat, Classes, K){

  require(dplyr)
  require(caret)
  require(glmnet)

  df <- data.frame(ID = colnames(Mat), Classes = Classes)
  samples <- createFolds(df$Classes, k = K, list = TRUE)
  ##stratified k-fold cross-validation
  #samples <- createFolds(df$Classes, k = K, list = TRUE, stratified = TRUE)
  ##Index in samples are test set indices
  return(list(df = df, samples = samples))

}

# The following script is only used for k-Fold cross-validation
OnevsEach.HCC <- function(Mat, classes.df, Indices, nDMR, subclass) {

  require(randomForest)

  # Training/test-set split

  TestData <- Mat[,Indices]
  TestPheno <- classes.df[Indices,]

  TrainData <- Mat[,!colnames(Mat) %in% TestPheno$ID]
  TrainPheno <- classes.df %>% filter(!ID %in% TestPheno$ID)

  DMRList <- list()

  # This loop does DMR preselection using a one vs each criterion

  FixedClass <- which(TrainPheno$Classes == subclass)
  OtherClasses <- which(TrainPheno$Classes != subclass)

  FixedClass.matrix <- TrainData[,FixedClass]
  OtherMatrix <- TrainData[,OtherClasses]
  DMR.classes <- c(rep(subclass,ncol(FixedClass.matrix)), rep("Others",ncol(OtherMatrix)))

  DMR.Data <- cbind(FixedClass.matrix, OtherMatrix)
  Des <- model.matrix(~0 + DMR.classes)
  colnames(Des) <- levels(factor(DMR.classes))

  LimmaFit <- lmFit(DMR.Data, Des) %>%
    contrasts.fit(., makeContrasts(eval(parse(text = paste0(subclass, "-Others"))), levels = Des)) %>%
    eBayes(., trend = TRUE) %>%
    topTable(., number = nrow(FixedClass.matrix))

  LimmaFit <- LimmaFit %>% .[order(.$t),]

  nDMR.b <- nDMR / 2

  TotalRows <- nrow(LimmaFit) - (nDMR.b - 1)
  Features <- rbind(LimmaFit[1:nDMR.b,], LimmaFit[TotalRows:nrow(LimmaFit),])

  Features <- rownames(Features)
  DMRList[[subclass]] <- Features
  message(length(Features))
  message("FixedClass:")
  message(paste0(length(FixedClass.matrix)," of FixedClass"))
  message(colnames(FixedClass.matrix))
  message("OtherMatrix:")
  message(paste0(length(OtherMatrix)," of OtherClass"))
  message(colnames(OtherMatrix))
  message(paste0(length(Features), " of DMR"))
  message(paste0(subclass," vs other classes DMR selection done"))

  # This creates feature set
  Features <- unlist(DMRList)

  # Here we fit the model and chuck it into the modlist
  NewAnn <- ifelse(TrainPheno$Classes == subclass,"One","Others")
  mtry.val <- nrow(TrainData[rownames(TrainData) %in% Features,]) / 3

  Model <- train(x = t(TrainData[rownames(TrainData) %in% Features,]), y = factor(NewAnn), trControl = Features.CVparam, method = "rf" , tuneGrid = expand.grid(.mtry = mtry.val), metric = "Kappa")
  message("Model Selection Complete")
  Prediction.classProbs <- predict(Model, newdata = t(TestData), type = "prob") %>%
    data.frame

  Prediction.classProbs$ActualClass <- TestPheno$Classes
  Prediction.classProbs$PredictedClass <- predict(Model, newdata = t(TestData), type = "raw")


  CombinedOutput <- list(Model = Model, TestPred = Prediction.classProbs)

  return(CombinedOutput)

}
