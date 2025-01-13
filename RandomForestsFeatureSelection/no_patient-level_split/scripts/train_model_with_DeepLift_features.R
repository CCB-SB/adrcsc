suppressMessages(library(Seurat))
suppressMessages(library(yaml))
suppressMessages(library(caret))
suppressMessages(library(ranger))
suppressMessages(library(magrittr))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(doParallel))
suppressMessages(library(ggplotify))

set.seed(123)

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile, type = "output")
sink(logfile, type = "message")

### Custom data class
DataClass <- setClass("DataClass", slots = c(x_train = "data.frame", y_train = "factor", x_v2 = "data.frame", y_v2 = "factor", x_test = "data.frame", y_test = "factor"))

### Helper functions
largest_char_vector <- function(...) {
  vectors <- list(...)
  return(vectors[[which.max(sapply(vectors, length))]])
}

### Load DeepLift feature ranking and define feature sets
DL_ranking <- read.csv(file = snakemake@input[[1]])
top30_features <- DL_ranking$X[1:30]
scGR_selected_features_defaultconfig <- read.csv(file = snakemake@input[[5]])
DL_features_scGR_defaultconfig <- DL_ranking$X[1:length(scGR_selected_features_defaultconfig$gene)]
scGR_selected_features_bestconfig <- read.csv(file = snakemake@input[[6]])
DL_features_scGR_bestconfig <- DL_ranking$X[1:length(scGR_selected_features_bestconfig$gene)]
largest_feature_set <- largest_char_vector(top30_features, DL_features_scGR_defaultconfig, DL_features_scGR_bestconfig)

### Load and subset data
load(snakemake@input[[2]])
seu <- readRDS(snakemake@input[[3]])
seu_train <- subset(seu, cells = c(rownames(data@x_train), rownames(data@x_v2)))
DefaultAssay(seu_train) <- "RNA"
seu_train <- NormalizeData(seu_train)

### Get expression values and labels
mtx <- seu_train@assays$RNA@data[largest_feature_set, ]
mtx <- t(mtx)
data@x_train <- as.data.frame(mtx)

common_levels <- union(levels(data@y_train), levels(data@y_v2))
y_trainV1V2 <- factor(c(as.character(data@y_train), as.character(data@y_v2)), levels = common_levels)
data@y_train <- y_trainV1V2

gc()

### Train RF model on full TrainV1V2 with different DeepLift selections
top30_features <- gsub("-", ".", top30_features)
DL_features_scGR_defaultconfig <- gsub("-", ".", DL_features_scGR_defaultconfig)
DL_features_scGR_bestconfig <- gsub("-", ".", DL_features_scGR_bestconfig)

xy_train <- cbind(data@x_train, y = data@y_train)
colnames(xy_train) <- gsub("-", ".", colnames(xy_train))

rf_model_top30 <- ranger(y ~ ., data = xy_train[, c(top30_features, "y")], importance = "impurity")
rf_model_scGR_defaultconfig <- ranger(y ~ ., data = xy_train[, c(DL_features_scGR_defaultconfig, "y")], importance = "impurity")
rf_model_scGR_bestconfig <- ranger(y ~ ., data = xy_train[, c(DL_features_scGR_bestconfig, "y")], importance = "impurity")

### Assess performance on test sets Test1 and Test2
seu_test1 <- subset(seu, cells = rownames(data@x_test))
seu_test2 <- readRDS(file = snakemake@input[[4]])
DefaultAssay(seu_test1) <- "RNA"
DefaultAssay(seu_test2) <- "RNA"
seu_test1 <- NormalizeData(seu_test1)
seu_test2 <- NormalizeData(seu_test2)

# Test1 (Test cells from same patients)
genes_with_dot <- setdiff(gsub(".", "-", largest_feature_set, fixed = TRUE), rownames(seu_test1@assays$RNA@data))
mtx_genes <- intersect(gsub(".", "-", largest_feature_set, fixed = TRUE), rownames(seu_test1@assays$RNA@data))
mtx_genes <- c(mtx_genes, gsub("-", ".", genes_with_dot))
mtx <- seu_test1@assays$RNA@data[mtx_genes, ]
mtx <- t(mtx)
data@x_test <- as.data.frame(mtx)
colnames(data@x_test) <- gsub("-", ".", colnames(data@x_test))
gc()

predictions <- predict(rf_model_top30, data = data@x_test[, top30_features], type = "response")$predictions
cat("Confusion matrix of top 30 DeepLift features RF model on Test1:", "\n")
cm <- confusionMatrix(predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_RF_Test1_top30 <- cm$overall[["Accuracy"]]
F1_RF_Test1_top30 <- cm$byClass[["F1"]]

predictions <- predict(rf_model_scGR_defaultconfig, data = data@x_test[, DL_features_scGR_defaultconfig], type = "response")$predictions
cat("Confusion matrix of top", length(DL_features_scGR_defaultconfig), "(set size from scGR defaultconfig) DeepLift features RF model on Test1:", "\n")
cm <- confusionMatrix(predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_RF_Test1_scGR_defaultconfig <- cm$overall[["Accuracy"]]
F1_RF_Test1_scGR_defaultconfig <- cm$byClass[["F1"]]

predictions <- predict(rf_model_scGR_bestconfig, data = data@x_test[, DL_features_scGR_bestconfig], type = "response")$predictions
cat("Confusion matrix of top", length(DL_features_scGR_bestconfig), "(set size from scGR bestconfig) DeepLift features RF model on Test1:", "\n")
cm <- confusionMatrix(predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_RF_Test1_scGR_bestconfig <- cm$overall[["Accuracy"]]
F1_RF_Test1_scGR_bestconfig <- cm$byClass[["F1"]]

# Test2 (Test cells from independent patients)
annotation_slot <- snakemake@config[["annotation_slot"]]
expr <- FetchData(object = seu_test2, vars = annotation_slot)
seu_test2 <- seu_test2[, which(x = expr == snakemake@config[["celltypes"]][[snakemake@params[[1]]]][["name"]])]
message(snakemake@params[[1]], " Test2 subset object: ")
seu_test2
message("# ", snakemake@config[["label_of_interest"]], " cells:")
message(sum(seu_test2[[snakemake@config[["class_of_interest"]]]] == snakemake@config[["label_of_interest"]]))
message("# other cells:")
message(sum(seu_test2[[snakemake@config[["class_of_interest"]]]] != snakemake@config[["label_of_interest"]]))

genes_with_dot <- setdiff(gsub(".", "-", largest_feature_set, fixed = TRUE), rownames(seu_test2@assays$RNA@data))
mtx_genes <- intersect(gsub(".", "-", largest_feature_set, fixed = TRUE), rownames(seu_test2@assays$RNA@data))
mtx_genes <- c(mtx_genes, gsub("-", ".", genes_with_dot))
mtx <- seu_test2@assays$RNA@data[mtx_genes, ]
mtx <- t(mtx)
x_test2 <- as.data.frame(mtx)
colnames(x_test2) <- gsub("-", ".", colnames(x_test2))
y_test2 <- seu_test2$Diagnosis
y_test2 <- ifelse(y_test2 == snakemake@config[["label_of_interest"]], snakemake@config[["label_of_interest"]], "ZZZ")
y_test2 <- factor(y_test2)

predictions <- predict(rf_model_top30, data = x_test2[, top30_features], type = "response")$predictions
cat("Confusion matrix of top 30 DeepLift features RF model on Test2:", "\n")
cm <- confusionMatrix(predictions, y_test2, mode = "everything")  
print(cm)
Accuracy_RF_Test2_top30 <- cm$overall[["Accuracy"]]
F1_RF_Test2_top30 <- cm$byClass[["F1"]]

predictions <- predict(rf_model_scGR_defaultconfig, data = x_test2[, DL_features_scGR_defaultconfig], type = "response")$predictions
cat("Confusion matrix of top", length(DL_features_scGR_defaultconfig), "(set size from scGR defaultconfig) DeepLift features RF model on Test2:", "\n")
cm <- confusionMatrix(predictions, y_test2, mode = "everything")  
print(cm)
Accuracy_RF_Test2_scGR_defaultconfig <- cm$overall[["Accuracy"]]
F1_RF_Test2_scGR_defaultconfig <- cm$byClass[["F1"]]

predictions <- predict(rf_model_scGR_bestconfig, data = x_test2[, DL_features_scGR_bestconfig], type = "response")$predictions
cat("Confusion matrix of top", length(DL_features_scGR_bestconfig), "(set size from scGR bestconfig) DeepLift features RF model on Test2:", "\n")
cm <- confusionMatrix(predictions, y_test2, mode = "everything")  
print(cm)
Accuracy_RF_Test2_scGR_bestconfig <- cm$overall[["Accuracy"]]
F1_RF_Test2_scGR_bestconfig <- cm$byClass[["F1"]]

### Train kNN classifier for comparison across feature selection methods
set.seed(123)
train_control <- trainControl(method = "cv", number = 10)
cl <- makePSOCKcluster(snakemake@threads)
registerDoParallel(cl)

start_time <- Sys.time()
knn_model_top30 <- train(y ~ ., data = xy_train[, c(top30_features, "y")], 
                         method = "knn", trControl = train_control,
                         tuneLength = 10)
end_time <- Sys.time()
print("kNN training (top30) finished after:")
print(end_time - start_time)
print(knn_model_top30)

start_time <- Sys.time()
knn_model_scGR_defaultconfig <- train(y ~ ., data = xy_train[, c(DL_features_scGR_defaultconfig, "y")], 
                                      method = "knn", trControl = train_control,
                                      tuneLength = 10)
end_time <- Sys.time()
print("kNN training (scGR defaultconfig set size) finished after:")
print(end_time - start_time)
print(knn_model_scGR_defaultconfig)

start_time <- Sys.time()
knn_model_scGR_bestconfig <- train(y ~ ., data = xy_train[, c(DL_features_scGR_bestconfig, "y")], 
                                   method = "knn", trControl = train_control,
                                   tuneLength = 10)
end_time <- Sys.time()
print("kNN training (scGR bestconfig set size) finished after:")
print(end_time - start_time)
print(knn_model_scGR_bestconfig)

stopCluster(cl)

knn_predictions <- predict(knn_model_top30, newdata = data@x_test[, top30_features])
cat("Confusion matrix of top 30 DeepLift features kNN model on Test1:", "\n")
cm <- confusionMatrix(knn_predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_kNN_Test1_top30 <- cm$overall[["Accuracy"]]
F1_kNN_Test1_top30 <- cm$byClass[["F1"]]
knn_predictions <- predict(knn_model_top30, newdata = x_test2[, top30_features])
cat("Confusion matrix of top 30 DeepLift features kNN model on Test2:", "\n")
cm <- confusionMatrix(predictions, y_test2, mode = "everything")  
print(cm)
Accuracy_kNN_Test2_top30 <- cm$overall[["Accuracy"]]
F1_kNN_Test2_top30 <- cm$byClass[["F1"]]

knn_predictions <- predict(knn_model_scGR_defaultconfig, newdata = data@x_test[, DL_features_scGR_defaultconfig])
cat("Confusion matrix of top", length(DL_features_scGR_defaultconfig), "(set size from scGR defaultconfig) DeepLift features kNN model on Test1:", "\n")
cm <- confusionMatrix(knn_predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_kNN_Test1_scGR_defaultconfig <- cm$overall[["Accuracy"]]
F1_kNN_Test1_scGR_defaultconfig <- cm$byClass[["F1"]]
knn_predictions <- predict(knn_model_scGR_defaultconfig, newdata = x_test2[, DL_features_scGR_defaultconfig])
cat("Confusion matrix of top", length(DL_features_scGR_defaultconfig), "(set size from scGR defaultconfig) DeepLift features kNN model on Test2:", "\n")
cm <- confusionMatrix(knn_predictions, y_test2, mode = "everything")  
print(cm)
Accuracy_kNN_Test2_scGR_defaultconfig <- cm$overall[["Accuracy"]]
F1_kNN_Test2_scGR_defaultconfig <- cm$byClass[["F1"]]

knn_predictions <- predict(knn_model_scGR_bestconfig, newdata = data@x_test[, DL_features_scGR_bestconfig])
cat("Confusion matrix of top", length(DL_features_scGR_bestconfig), "(set size from scGR bestconfig) DeepLift features kNN model on Test1:", "\n")
cm <- confusionMatrix(knn_predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_kNN_Test1_scGR_bestconfig <- cm$overall[["Accuracy"]]
F1_kNN_Test1_scGR_bestconfig <- cm$byClass[["F1"]]
knn_predictions <- predict(knn_model_scGR_bestconfig, newdata = x_test2[, DL_features_scGR_bestconfig])
cat("Confusion matrix of top", length(DL_features_scGR_bestconfig), "(set size from scGR bestconfig) DeepLift features kNN model on Test2:", "\n")
cm <- confusionMatrix(knn_predictions, y_test2, mode = "everything")  
print(cm)
Accuracy_kNN_Test2_scGR_bestconfig <- cm$overall[["Accuracy"]]
F1_kNN_Test2_scGR_bestconfig <- cm$byClass[["F1"]]

### Patient-level predictions
donor_group <- c()
donor_slot <- "SCMD"
meta_data_column <- "Diagnosis"
group2 <- snakemake@config[["label_of_interest"]]
donors <- unique(seu_test2[[donor_slot]][, 1])
for (donor in donors) {
  if (sum(seu_test2[[donor_slot]] == donor & seu_test2[[meta_data_column]] == group2) > 0) {
    donor_group <- append(donor_group, 1)
  }
  else {
    donor_group <- append(donor_group, 0)
  }
}
donor_df <- as.data.frame(cbind(donors, donor_group))
colnames(donor_df) <- c("Donor", "Group.TrueLabel")

donor_predictions_rf_top30 <- c()
prediction_percentages_rf_top30 <- c()
donor_predictions_knn_top30 <- c()
prediction_percentages_knn_top30 <- c()

donor_predictions_rf_scGR_defaultconfig <- c()
prediction_percentages_rf_scGR_defaultconfig <- c()
donor_predictions_knn_scGR_defaultconfig <- c()
prediction_percentages_knn_scGR_defaultconfig <- c()

donor_predictions_rf_scGR_bestconfig <- c()
prediction_percentages_rf_scGR_bestconfig <- c()
donor_predictions_knn_scGR_bestconfig <- c()
prediction_percentages_knn_scGR_bestconfig <- c()

for (donor in donors) {
  expr <- FetchData(object = seu_test2, vars = donor_slot)
  seu_donor <- seu_test2[, which(x = expr == donor)]
  if (ncol(seu_donor) < 2) {
    donor_df <- donor_df[donor_df$Donor != donor, ]
    next
  }

  mtx <- seu_donor@assays$RNA@data[mtx_genes, ]
  mtx <- t(mtx)
  x_donor <- as.data.frame(mtx)
  colnames(x_donor) <- gsub("-", ".", colnames(x_donor))
  y_donor <- seu_donor$Diagnosis
  y_donor <- ifelse(y_donor == snakemake@config[["label_of_interest"]], snakemake@config[["label_of_interest"]], "ZZZ")
  y_donor <- factor(y_donor)

  rf_pred <- predict(rf_model_top30, data = x_donor[, top30_features], type = "response")$predictions
  pos_pred <- sum(rf_pred == group2)
  neg_pred <- sum(rf_pred == "ZZZ")
  if (pos_pred > neg_pred) {
    donor_predictions_rf_top30 <- append(donor_predictions_rf_top30, 1)
  }
  else {
    donor_predictions_rf_top30 <- append(donor_predictions_rf_top30, 0)
  }
  prediction_percentages_rf_top30 <- append(prediction_percentages_rf_top30, pos_pred / (pos_pred + neg_pred))
  knn_pred <- predict(knn_model_top30, newdata = x_donor[, top30_features])
  pos_pred <- sum(knn_pred == group2)
  neg_pred <- sum(knn_pred == "ZZZ")
  if (pos_pred > neg_pred) {
    donor_predictions_knn_top30 <- append(donor_predictions_knn_top30, 1)
  }
  else {
    donor_predictions_knn_top30 <- append(donor_predictions_knn_top30, 0)
  }
  prediction_percentages_knn_top30 <- append(prediction_percentages_knn_top30, pos_pred / (pos_pred + neg_pred))

  rf_pred <- predict(rf_model_scGR_defaultconfig, data = x_donor[, DL_features_scGR_defaultconfig], type = "response")$predictions
  pos_pred <- sum(rf_pred == group2)
  neg_pred <- sum(rf_pred == "ZZZ")
  if (pos_pred > neg_pred) {
    donor_predictions_rf_scGR_defaultconfig <- append(donor_predictions_rf_scGR_defaultconfig, 1)
  }
  else {
    donor_predictions_rf_scGR_defaultconfig <- append(donor_predictions_rf_scGR_defaultconfig, 0)
  }
  prediction_percentages_rf_scGR_defaultconfig <- append(prediction_percentages_rf_scGR_defaultconfig, pos_pred / (pos_pred + neg_pred))
  knn_pred <- predict(knn_model_scGR_defaultconfig, newdata = x_donor[, DL_features_scGR_defaultconfig])
  pos_pred <- sum(knn_pred == group2)
  neg_pred <- sum(knn_pred == "ZZZ")
  if (pos_pred > neg_pred) {
    donor_predictions_knn_scGR_defaultconfig <- append(donor_predictions_knn_scGR_defaultconfig, 1)
  }
  else {
    donor_predictions_knn_scGR_defaultconfig <- append(donor_predictions_knn_scGR_defaultconfig, 0)
  }
  prediction_percentages_knn_scGR_defaultconfig <- append(prediction_percentages_knn_scGR_defaultconfig, pos_pred / (pos_pred + neg_pred))

  rf_pred <- predict(rf_model_scGR_bestconfig, data = x_donor[, DL_features_scGR_bestconfig], type = "response")$predictions
  pos_pred <- sum(rf_pred == group2)
  neg_pred <- sum(rf_pred == "ZZZ")
  if (pos_pred > neg_pred) {
    donor_predictions_rf_scGR_bestconfig <- append(donor_predictions_rf_scGR_bestconfig, 1)
  }
  else {
    donor_predictions_rf_scGR_bestconfig <- append(donor_predictions_rf_scGR_bestconfig, 0)
  }
  prediction_percentages_rf_scGR_bestconfig <- append(prediction_percentages_rf_scGR_bestconfig, pos_pred / (pos_pred + neg_pred))
  knn_pred <- predict(knn_model_scGR_bestconfig, newdata = x_donor[, DL_features_scGR_bestconfig])
  pos_pred <- sum(knn_pred == group2)
  neg_pred <- sum(knn_pred == "ZZZ")
  if (pos_pred > neg_pred) {
    donor_predictions_knn_scGR_bestconfig <- append(donor_predictions_knn_scGR_bestconfig, 1)
  }
  else {
    donor_predictions_knn_scGR_bestconfig <- append(donor_predictions_knn_scGR_bestconfig, 0)
  }
  prediction_percentages_knn_scGR_bestconfig <- append(prediction_percentages_knn_scGR_bestconfig, pos_pred / (pos_pred + neg_pred))
}

donor_df <- cbind(donor_df, donor_predictions_rf_top30, prediction_percentages_rf_top30, donor_predictions_knn_top30, prediction_percentages_knn_top30)
donor_df <- cbind(donor_df, donor_predictions_rf_scGR_defaultconfig, prediction_percentages_rf_scGR_defaultconfig, donor_predictions_knn_scGR_defaultconfig, prediction_percentages_knn_scGR_defaultconfig)
donor_df <- cbind(donor_df, donor_predictions_rf_scGR_bestconfig, prediction_percentages_rf_scGR_bestconfig, donor_predictions_knn_scGR_bestconfig, prediction_percentages_knn_scGR_bestconfig)
write.csv(donor_df, file = file.path(snakemake@output[[2]]), row.names = FALSE)

cat("Patient-level predictions: ", "\n")

cm_donor_pred <- confusionMatrix(factor(donor_predictions_rf_top30, levels = c(1, 0)), factor(donor_df$Group.TrueLabel, levels = c(1, 0)), mode = "everything")
cat("Confusion matrix of top 30 DeepLift features RF model patient-level predictions on Test2:", "\n")
print(cm_donor_pred)
Accuracy_RF_Test2_Donor_Prediction_top30 <- cm_donor_pred$overall[["Accuracy"]]
F1_RF_Test2_Donor_Prediction_top30 <- cm_donor_pred$byClass[["F1"]]
cm_donor_pred <- confusionMatrix(factor(donor_predictions_knn_top30, levels = c(1, 0)), factor(donor_df$Group.TrueLabel, levels = c(1, 0)), mode = "everything")
cat("Confusion matrix of top 30 DeepLift features kNN model patient-level predictions on Test2:", "\n")
print(cm_donor_pred)
Accuracy_kNN_Test2_Donor_Prediction_top30 <- cm_donor_pred$overall[["Accuracy"]]
F1_kNN_Test2_Donor_Prediction_top30 <- cm_donor_pred$byClass[["F1"]]

cm_donor_pred <- confusionMatrix(factor(donor_predictions_rf_scGR_defaultconfig, levels = c(1, 0)), factor(donor_df$Group.TrueLabel, levels = c(1, 0)), mode = "everything")
cat("Confusion matrix of top", length(DL_features_scGR_defaultconfig), "(set size from scGR defaultconfig) DeepLift features RF model patient-level predictions on Test2:", "\n")
print(cm_donor_pred)
Accuracy_RF_Test2_Donor_Prediction_scGR_defaultconfig <- cm_donor_pred$overall[["Accuracy"]]
F1_RF_Test2_Donor_Prediction_scGR_defaultconfig <- cm_donor_pred$byClass[["F1"]]
cm_donor_pred <- confusionMatrix(factor(donor_predictions_knn_scGR_defaultconfig, levels = c(1, 0)), factor(donor_df$Group.TrueLabel, levels = c(1, 0)), mode = "everything")
cat("Confusion matrix of top", length(DL_features_scGR_defaultconfig), "(set size from scGR defaultconfig) DeepLift features kNN model patient-level predictions on Test2:", "\n")
print(cm_donor_pred)
Accuracy_kNN_Test2_Donor_Prediction_scGR_defaultconfig <- cm_donor_pred$overall[["Accuracy"]]
F1_kNN_Test2_Donor_Prediction_scGR_defaultconfig <- cm_donor_pred$byClass[["F1"]]

cm_donor_pred <- confusionMatrix(factor(donor_predictions_rf_scGR_bestconfig, levels = c(1, 0)), factor(donor_df$Group.TrueLabel, levels = c(1, 0)), mode = "everything")
cat("Confusion matrix of top", length(DL_features_scGR_bestconfig), "(set size from scGR bestconfig) DeepLift features RF model patient-level predictions on Test2:", "\n")
print(cm_donor_pred)
Accuracy_RF_Test2_Donor_Prediction_scGR_bestconfig <- cm_donor_pred$overall[["Accuracy"]]
F1_RF_Test2_Donor_Prediction_scGR_bestconfig <- cm_donor_pred$byClass[["F1"]]
cm_donor_pred <- confusionMatrix(factor(donor_predictions_knn_scGR_bestconfig, levels = c(1, 0)), factor(donor_df$Group.TrueLabel, levels = c(1, 0)), mode = "everything")
cat("Confusion matrix of top", length(DL_features_scGR_bestconfig), "(set size from scGR bestconfig) DeepLift features kNN model patient-level predictions on Test2:", "\n")
print(cm_donor_pred)
Accuracy_kNN_Test2_Donor_Prediction_scGR_bestconfig <- cm_donor_pred$overall[["Accuracy"]]
F1_kNN_Test2_Donor_Prediction_scGR_bestconfig <- cm_donor_pred$byClass[["F1"]]

### Save results
DL_ranking <- DL_ranking[1:length(largest_feature_set), ]

DL_ranking$Accuracy_RF_Test1_top30 <- Accuracy_RF_Test1_top30
DL_ranking$Accuracy_RF_Test2_top30 <- Accuracy_RF_Test2_top30
DL_ranking$Accuracy_kNN_Test1_top30 <- Accuracy_kNN_Test1_top30
DL_ranking$Accuracy_kNN_Test2_top30 <- Accuracy_kNN_Test2_top30
DL_ranking$F1_RF_Test1_top30 <- F1_RF_Test1_top30
DL_ranking$F1_RF_Test2_top30 <- F1_RF_Test2_top30
DL_ranking$F1_kNN_Test1_top30 <- F1_kNN_Test1_top30
DL_ranking$F1_kNN_Test2_top30 <- F1_kNN_Test2_top30
DL_ranking$Accuracy_RF_Test2_Donor_Prediction_top30 <- Accuracy_RF_Test2_Donor_Prediction_top30
DL_ranking$Accuracy_kNN_Test2_Donor_Prediction_top30 <- Accuracy_kNN_Test2_Donor_Prediction_top30
DL_ranking$F1_RF_Test2_Donor_Prediction_top30 <- F1_RF_Test2_Donor_Prediction_top30
DL_ranking$F1_kNN_Test2_Donor_Prediction_top30 <- F1_kNN_Test2_Donor_Prediction_top30

DL_ranking$Accuracy_RF_Test1_scGR_defaultconfig <- Accuracy_RF_Test1_scGR_defaultconfig
DL_ranking$Accuracy_RF_Test2_scGR_defaultconfig <- Accuracy_RF_Test2_scGR_defaultconfig
DL_ranking$Accuracy_kNN_Test1_scGR_defaultconfig <- Accuracy_kNN_Test1_scGR_defaultconfig
DL_ranking$Accuracy_kNN_Test2_scGR_defaultconfig <- Accuracy_kNN_Test2_scGR_defaultconfig
DL_ranking$F1_RF_Test1_scGR_defaultconfig <- F1_RF_Test1_scGR_defaultconfig
DL_ranking$F1_RF_Test2_scGR_defaultconfig <- F1_RF_Test2_scGR_defaultconfig
DL_ranking$F1_kNN_Test1_scGR_defaultconfig <- F1_kNN_Test1_scGR_defaultconfig
DL_ranking$F1_kNN_Test2_scGR_defaultconfig <- F1_kNN_Test2_scGR_defaultconfig
DL_ranking$Accuracy_RF_Test2_Donor_Prediction_scGR_defaultconfig <- Accuracy_RF_Test2_Donor_Prediction_scGR_defaultconfig
DL_ranking$Accuracy_kNN_Test2_Donor_Prediction_scGR_defaultconfig <- Accuracy_kNN_Test2_Donor_Prediction_scGR_defaultconfig
DL_ranking$F1_RF_Test2_Donor_Prediction_scGR_defaultconfig <- F1_RF_Test2_Donor_Prediction_scGR_defaultconfig
DL_ranking$F1_kNN_Test2_Donor_Prediction_scGR_defaultconfig <- F1_kNN_Test2_Donor_Prediction_scGR_defaultconfig

DL_ranking$Accuracy_RF_Test1_scGR_bestconfig <- Accuracy_RF_Test1_scGR_bestconfig
DL_ranking$Accuracy_RF_Test2_scGR_bestconfig <- Accuracy_RF_Test2_scGR_bestconfig
DL_ranking$Accuracy_kNN_Test1_scGR_bestconfig <- Accuracy_kNN_Test1_scGR_bestconfig
DL_ranking$Accuracy_kNN_Test2_scGR_bestconfig <- Accuracy_kNN_Test2_scGR_bestconfig
DL_ranking$F1_RF_Test1_scGR_bestconfig <- F1_RF_Test1_scGR_bestconfig
DL_ranking$F1_RF_Test2_scGR_bestconfig <- F1_RF_Test2_scGR_bestconfig
DL_ranking$F1_kNN_Test1_scGR_bestconfig <- F1_kNN_Test1_scGR_bestconfig
DL_ranking$F1_kNN_Test2_scGR_bestconfig <- F1_kNN_Test2_scGR_bestconfig
DL_ranking$Accuracy_RF_Test2_Donor_Prediction_scGR_bestconfig <- Accuracy_RF_Test2_Donor_Prediction_scGR_bestconfig
DL_ranking$Accuracy_kNN_Test2_Donor_Prediction_scGR_bestconfig <- Accuracy_kNN_Test2_Donor_Prediction_scGR_bestconfig
DL_ranking$F1_RF_Test2_Donor_Prediction_scGR_bestconfig <- F1_RF_Test2_Donor_Prediction_scGR_bestconfig
DL_ranking$F1_kNN_Test2_Donor_Prediction_scGR_bestconfig <- F1_kNN_Test2_Donor_Prediction_scGR_bestconfig

write.csv(DL_ranking, file = snakemake@output[[1]], row.names = FALSE)