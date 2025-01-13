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
DataClass <- setClass("DataClass", slots = c(x_train = "data.frame", y_train = "factor", x_test = "data.frame", y_test = "factor"))

### Helper functions
largest_char_vector <- function(...) {
  vectors <- list(...)
  return(vectors[[which.max(sapply(vectors, length))]])
}

### Fetch values from config
class_column <- snakemake@config[["class_column"]]
annotation_column <- snakemake@config[["annotation_column"]]
donor_column <- snakemake@config[["donor_column"]]
positive_class <- snakemake@config[["positive_class"]]
celltype <- snakemake@config[["celltypes"]][[snakemake@params[[1]]]][["name"]]

### Load DeepLift feature ranking and define feature sets
DL_ranking <- read.csv(file = snakemake@input[[1]])
top30_features <- DL_ranking$X[1:30]
scGR_selected_features_defaultconfig <- read.csv(file = snakemake@input[[5]])
DL_features_scGR_defaultconfig <- DL_ranking$X[1:length(scGR_selected_features_defaultconfig$gene)]
largest_feature_set <- largest_char_vector(top30_features, DL_features_scGR_defaultconfig)

### Load and subset data
load(snakemake@input[[2]])
seu <- readRDS(snakemake@input[[3]])
#seu_train <- subset(seu, cells = rownames(data@x_train))
expr <- FetchData(object = seu, vars = annotation_column)
seu_train <- seu[, which(x = expr == celltype)]
DefaultAssay(seu_train) <- "RNA"
seu_train <- NormalizeData(seu_train)

### Get expression values and labels
mtx <- seu_train@assays$RNA@data[largest_feature_set, ]
mtx <- t(mtx)
data@x_train <- as.data.frame(mtx)

y <- seu_train[[class_column]]
y <- ifelse(y == positive_class, positive_class, "ZZZ")
data@y_train <- factor(y)


gc()

### Train RF model on full TrainV1V2 with different DeepLift selections
top30_features <- gsub("-", ".", top30_features)
DL_features_scGR_defaultconfig <- gsub("-", ".", DL_features_scGR_defaultconfig)

set.seed(42)
xy_train <- cbind(data@x_train, y = data@y_train)
colnames(xy_train) <- gsub("-", ".", colnames(xy_train))
balanced_train_data <- downSample(xy_train, data@y_train, yname = "y")
dim(balanced_train_data)

rf_model_top30 <- ranger(y ~ ., data = balanced_train_data[, c(top30_features, "y")], importance = "impurity")
rf_model_scGR_defaultconfig <- ranger(y ~ ., data = balanced_train_data[, c(DL_features_scGR_defaultconfig, "y")], importance = "impurity")

### Assess performance on test sets Test1 and Test2
seu_test <- readRDS(file = snakemake@input[[4]])
seu_test <- subset(seu_test, cells = rownames(data@x_test))
DefaultAssay(seu_test) <- "RNA"
seu_test <- NormalizeData(seu_test)

genes_with_dot <- setdiff(gsub(".", "-", largest_feature_set, fixed = TRUE), rownames(seu_test@assays$RNA@data))
mtx_genes <- intersect(gsub(".", "-", largest_feature_set, fixed = TRUE), rownames(seu_test@assays$RNA@data))
mtx_genes <- c(mtx_genes, gsub("-", ".", genes_with_dot))
mtx <- seu_test@assays$RNA@data[mtx_genes, ]
mtx <- t(mtx)
data@x_test <- as.data.frame(mtx)
colnames(data@x_test) <- gsub("-", ".", colnames(data@x_test))
gc()

predictions <- predict(rf_model_top30, data = data@x_test[, top30_features], type = "response")$predictions
cat("Confusion matrix of top 30 DeepLift features RF model on Test set:", "\n")
cm <- confusionMatrix(predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_RF_top30 <- cm$overall[["Accuracy"]]
F1_RF_top30 <- cm$byClass[["F1"]]

predictions <- predict(rf_model_scGR_defaultconfig, data = data@x_test[, DL_features_scGR_defaultconfig], type = "response")$predictions
cat("Confusion matrix of top", length(DL_features_scGR_defaultconfig), "(set size from scGR defaultconfig) DeepLift features RF model on Test set:", "\n")
cm <- confusionMatrix(predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_RF_scGR_defaultconfig <- cm$overall[["Accuracy"]]
F1_RF_scGR_defaultconfig <- cm$byClass[["F1"]]

### Train kNN classifier for comparison across feature selection methods
set.seed(123)
train_control <- trainControl(method = "cv", number = 10)
cl <- makePSOCKcluster(snakemake@threads)
registerDoParallel(cl)

start_time <- Sys.time()
knn_model_top30 <- train(y ~ ., data = balanced_train_data[, c(top30_features, "y")], 
                         method = "knn", trControl = train_control,
                         tuneLength = 10)
end_time <- Sys.time()
print("kNN training (top30) finished after:")
print(end_time - start_time)
print(knn_model_top30)

start_time <- Sys.time()
knn_model_scGR_defaultconfig <- train(y ~ ., data = balanced_train_data[, c(DL_features_scGR_defaultconfig, "y")], 
                                      method = "knn", trControl = train_control,
                                      tuneLength = 10)
end_time <- Sys.time()
print("kNN training (scGR defaultconfig set size) finished after:")
print(end_time - start_time)
print(knn_model_scGR_defaultconfig)

stopCluster(cl)

knn_predictions <- predict(knn_model_top30, newdata = data@x_test[, top30_features])
cat("Confusion matrix of top 30 DeepLift features kNN model on Test set:", "\n")
cm <- confusionMatrix(knn_predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_kNN_top30 <- cm$overall[["Accuracy"]]
F1_kNN_top30 <- cm$byClass[["F1"]]

knn_predictions <- predict(knn_model_scGR_defaultconfig, newdata = data@x_test[, DL_features_scGR_defaultconfig])
cat("Confusion matrix of top", length(DL_features_scGR_defaultconfig), "(set size from scGR defaultconfig) DeepLift features kNN model on Test set:", "\n")
cm <- confusionMatrix(knn_predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_kNN_scGR_defaultconfig <- cm$overall[["Accuracy"]]
F1_kNN_scGR_defaultconfig <- cm$byClass[["F1"]]

### Patient-level predictions
donor_group <- c()
donors <- unique(seu_test[[donor_column]][, 1])
for (donor in donors) {
  if (sum(seu_test[[donor_column]] == donor & seu_test[[class_column]] == positive_class) > 0) {
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

for (donor in donors) {
  expr <- FetchData(object = seu_test, vars = donor_column)
  seu_donor <- seu_test[, which(x = expr == donor)]
  if (ncol(seu_donor) < 2) {
    donor_df <- donor_df[donor_df$Donor != donor, ]
    next
  }

  mtx <- seu_donor@assays$RNA@data[mtx_genes, ]
  mtx <- t(mtx)
  x_donor <- as.data.frame(mtx)
  colnames(x_donor) <- gsub("-", ".", colnames(x_donor))
  y_donor <- seu_donor[[class_column]]
  y_donor <- ifelse(y_donor == positive_class, positive_class, "ZZZ")
  y_donor <- factor(y_donor)

  rf_pred <- predict(rf_model_top30, data = x_donor[, top30_features], type = "response")$predictions
  pos_pred <- sum(rf_pred == positive_class)
  neg_pred <- sum(rf_pred == "ZZZ")
  if (pos_pred > neg_pred) {
    donor_predictions_rf_top30 <- append(donor_predictions_rf_top30, 1)
  }
  else {
    donor_predictions_rf_top30 <- append(donor_predictions_rf_top30, 0)
  }
  prediction_percentages_rf_top30 <- append(prediction_percentages_rf_top30, pos_pred / (pos_pred + neg_pred))
  knn_pred <- predict(knn_model_top30, newdata = x_donor[, top30_features])
  pos_pred <- sum(knn_pred == positive_class)
  neg_pred <- sum(knn_pred == "ZZZ")
  if (pos_pred > neg_pred) {
    donor_predictions_knn_top30 <- append(donor_predictions_knn_top30, 1)
  }
  else {
    donor_predictions_knn_top30 <- append(donor_predictions_knn_top30, 0)
  }
  prediction_percentages_knn_top30 <- append(prediction_percentages_knn_top30, pos_pred / (pos_pred + neg_pred))

  rf_pred <- predict(rf_model_scGR_defaultconfig, data = x_donor[, DL_features_scGR_defaultconfig], type = "response")$predictions
  pos_pred <- sum(rf_pred == positive_class)
  neg_pred <- sum(rf_pred == "ZZZ")
  if (pos_pred > neg_pred) {
    donor_predictions_rf_scGR_defaultconfig <- append(donor_predictions_rf_scGR_defaultconfig, 1)
  }
  else {
    donor_predictions_rf_scGR_defaultconfig <- append(donor_predictions_rf_scGR_defaultconfig, 0)
  }
  prediction_percentages_rf_scGR_defaultconfig <- append(prediction_percentages_rf_scGR_defaultconfig, pos_pred / (pos_pred + neg_pred))
  knn_pred <- predict(knn_model_scGR_defaultconfig, newdata = x_donor[, DL_features_scGR_defaultconfig])
  pos_pred <- sum(knn_pred == positive_class)
  neg_pred <- sum(knn_pred == "ZZZ")
  if (pos_pred > neg_pred) {
    donor_predictions_knn_scGR_defaultconfig <- append(donor_predictions_knn_scGR_defaultconfig, 1)
  }
  else {
    donor_predictions_knn_scGR_defaultconfig <- append(donor_predictions_knn_scGR_defaultconfig, 0)
  }
  prediction_percentages_knn_scGR_defaultconfig <- append(prediction_percentages_knn_scGR_defaultconfig, pos_pred / (pos_pred + neg_pred))
}

donor_df <- cbind(donor_df, donor_predictions_rf_top30, prediction_percentages_rf_top30, donor_predictions_knn_top30, prediction_percentages_knn_top30)
donor_df <- cbind(donor_df, donor_predictions_rf_scGR_defaultconfig, prediction_percentages_rf_scGR_defaultconfig, donor_predictions_knn_scGR_defaultconfig, prediction_percentages_knn_scGR_defaultconfig)
write.csv(donor_df, file = file.path(snakemake@output[[2]]), row.names = FALSE)

cat("Patient-level predictions: ", "\n")

cm_donor_pred <- confusionMatrix(factor(donor_predictions_rf_top30, levels = c(1, 0)), factor(donor_df$Group.TrueLabel, levels = c(1, 0)), mode = "everything")
cat("Confusion matrix of top 30 DeepLift features RF model patient-level predictions on Test set", "\n")
print(cm_donor_pred)
Accuracy_RF_Test2_Donor_Prediction_top30 <- cm_donor_pred$overall[["Accuracy"]]
F1_RF_Test2_Donor_Prediction_top30 <- cm_donor_pred$byClass[["F1"]]
cm_donor_pred <- confusionMatrix(factor(donor_predictions_knn_top30, levels = c(1, 0)), factor(donor_df$Group.TrueLabel, levels = c(1, 0)), mode = "everything")
cat("Confusion matrix of top 30 DeepLift features kNN model patient-level predictions on Test set:", "\n")
print(cm_donor_pred)
Accuracy_kNN_Test2_Donor_Prediction_top30 <- cm_donor_pred$overall[["Accuracy"]]
F1_kNN_Test2_Donor_Prediction_top30 <- cm_donor_pred$byClass[["F1"]]

cm_donor_pred <- confusionMatrix(factor(donor_predictions_rf_scGR_defaultconfig, levels = c(1, 0)), factor(donor_df$Group.TrueLabel, levels = c(1, 0)), mode = "everything")
cat("Confusion matrix of top", length(DL_features_scGR_defaultconfig), "(set size from scGR defaultconfig) DeepLift features RF model patient-level predictions on Test set:", "\n")
print(cm_donor_pred)
Accuracy_RF_Test2_Donor_Prediction_scGR_defaultconfig <- cm_donor_pred$overall[["Accuracy"]]
F1_RF_Test2_Donor_Prediction_scGR_defaultconfig <- cm_donor_pred$byClass[["F1"]]
cm_donor_pred <- confusionMatrix(factor(donor_predictions_knn_scGR_defaultconfig, levels = c(1, 0)), factor(donor_df$Group.TrueLabel, levels = c(1, 0)), mode = "everything")
cat("Confusion matrix of top", length(DL_features_scGR_defaultconfig), "(set size from scGR defaultconfig) DeepLift features kNN model patient-level predictions on Test set:", "\n")
print(cm_donor_pred)
Accuracy_kNN_Test2_Donor_Prediction_scGR_defaultconfig <- cm_donor_pred$overall[["Accuracy"]]
F1_kNN_Test2_Donor_Prediction_scGR_defaultconfig <- cm_donor_pred$byClass[["F1"]]

### Save results
DL_ranking <- DL_ranking[1:length(largest_feature_set), ]

DL_ranking$Accuracy_RF_top30 <- Accuracy_RF_top30
DL_ranking$Accuracy_kNN_top30 <- Accuracy_kNN_top30
DL_ranking$F1_RF_top30 <- F1_RF_top30 
DL_ranking$F1_kNN_top30 <- F1_kNN_top30
DL_ranking$Accuracy_RF_Test2_Donor_Prediction_top30 <- Accuracy_RF_Test2_Donor_Prediction_top30
DL_ranking$Accuracy_kNN_Test2_Donor_Prediction_top30 <- Accuracy_kNN_Test2_Donor_Prediction_top30
DL_ranking$F1_RF_Test2_Donor_Prediction_top30 <- F1_RF_Test2_Donor_Prediction_top30
DL_ranking$F1_kNN_Test2_Donor_Prediction_top30 <- F1_kNN_Test2_Donor_Prediction_top30

DL_ranking$Accuracy_RF_scGR_defaultconfig <- Accuracy_RF_scGR_defaultconfig
DL_ranking$Accuracy_kNN_scGR_defaultconfig <- Accuracy_kNN_scGR_defaultconfig
DL_ranking$F1_RF_scGR_defaultconfig <- F1_RF_scGR_defaultconfig
DL_ranking$F1_kNN_scGR_defaultconfig <- F1_kNN_scGR_defaultconfig
DL_ranking$Accuracy_RF_Test2_Donor_Prediction_scGR_defaultconfig <- Accuracy_RF_Test2_Donor_Prediction_scGR_defaultconfig
DL_ranking$Accuracy_kNN_Test2_Donor_Prediction_scGR_defaultconfig <- Accuracy_kNN_Test2_Donor_Prediction_scGR_defaultconfig
DL_ranking$F1_RF_Test2_Donor_Prediction_scGR_defaultconfig <- F1_RF_Test2_Donor_Prediction_scGR_defaultconfig
DL_ranking$F1_kNN_Test2_Donor_Prediction_scGR_defaultconfig <- F1_kNN_Test2_Donor_Prediction_scGR_defaultconfig

write.csv(DL_ranking, file = snakemake@output[[1]], row.names = FALSE)