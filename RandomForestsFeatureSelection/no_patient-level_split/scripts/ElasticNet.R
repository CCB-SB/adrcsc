suppressMessages(library(Seurat))
suppressMessages(library(yaml))
suppressMessages(library(caret))
suppressMessages(library(ranger))
suppressMessages(library(glmnet))
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

### scElasticNetFS
### Function to perform feature selection via elastic net logistic regression
sc_ElasticNet_FS <- function(x, y, alpha = 0.5, nfolds = 5) {
  
  x <- as.matrix(x)
  y <- ifelse(y == "ZZZ", 0, 1) #factor?
  
  cl <- makePSOCKcluster(snakemake@threads)
  registerDoParallel(cl)

  start_time <- Sys.time()
  cvfit <- cv.glmnet(x, y, family = "binomial", alpha = alpha, type.measure = "class", nfolds = nfolds, standardize = TRUE, parallel=TRUE)
  end_time <- Sys.time()

  stopCluster(cl)
  
  elasticnet_time <- end_time - start_time
  print(elasticnet_time)
  
  cvfit
}

### Load and subset data
load(snakemake@input[[1]])
seu <- readRDS(snakemake@input[[2]])
seu_train <- subset(seu, cells = c(rownames(data@x_train), rownames(data@x_v2)))
DefaultAssay(seu_train) <- "RNA"
seu_train <- NormalizeData(seu_train)

### Filtering of lowly expressed genes
expr_matrix <- GetAssayData(seu_train, slot = "data")
seu_train <- SetIdent(seu_train, value = "Diagnosis")
conditions <- unique(seu_train$Diagnosis)

# Initialize a logical vector to track genes that are expressed in >=1% of cells in at least one condition
genes_to_keep <- rep(FALSE, nrow(expr_matrix))
names(genes_to_keep) <- rownames(expr_matrix)

# Loop over each condition and apply the filter
for (condition in conditions) {
  
  # Subset cells belonging to the current condition
  condition_cells <- WhichCells(seu_train, ident = condition)
  
  # Subset the expression matrix for the current condition
  expr_data_condition <- expr_matrix[, condition_cells]
  
  # Calculate the percentage of cells in which each gene is expressed
  percent_expressed <- Matrix::rowSums(expr_data_condition > 0) / length(condition_cells)
  
  # Keep genes that are expressed in >=1% of cells in this condition
  genes_to_keep <- genes_to_keep | (percent_expressed >= 0.01)
}

input_genes <- names(genes_to_keep[genes_to_keep])


### Get expression values and labels
mtx <- seu_train@assays$RNA@data[input_genes, ]
mtx <- t(mtx)
data@x_train <- as.data.frame(mtx)

common_levels <- union(levels(data@y_train), levels(data@y_v2))
y_trainV1V2 <- factor(c(as.character(data@y_train), as.character(data@y_v2)), levels = common_levels)
data@y_train <- y_trainV1V2

gc()

### Run Logistic Regression with ElasticNet Regularization (alpha = 0.5)
xy_train <- cbind(data@x_train, y = data@y_train)
colnames(xy_train) <- gsub("-", ".", colnames(xy_train))
colnames(data@x_train) <- gsub("-", ".", colnames(data@x_train))

sc_ElasticNet_FS_results <- sc_ElasticNet_FS(data@x_train, data@y_train)
feature_ranking_ElasticNet <- coef(sc_ElasticNet_FS_results, s = "lambda.1se")
feature_ranking_ElasticNet <- feature_ranking_ElasticNet[order(feature_ranking_ElasticNet, decreasing = TRUE), ]
selected_features_ElasticNet <- feature_ranking_ElasticNet[feature_ranking_ElasticNet != 0]
selected_features_ElasticNet <- selected_features_ElasticNet[names(selected_features_ElasticNet) != "(Intercept)"]

### Assess performance on test sets Test1 and Test2
seu_test1 <- subset(seu, cells = rownames(data@x_test))
seu_test2 <- readRDS(file = snakemake@input[[3]])
DefaultAssay(seu_test1) <- "RNA"
DefaultAssay(seu_test2) <- "RNA"
seu_test1 <- NormalizeData(seu_test1)
seu_test2 <- NormalizeData(seu_test2)

# Test1 (Test cells from same patients)
genes_with_dot <- setdiff(gsub(".", "-", input_genes, fixed = TRUE), rownames(seu_test1@assays$RNA@data))
mtx_genes <- intersect(gsub(".", "-", input_genes, fixed = TRUE), rownames(seu_test1@assays$RNA@data))
mtx_genes <- c(mtx_genes, gsub("-", ".", genes_with_dot))
mtx <- seu_test1@assays$RNA@data[mtx_genes, ]
mtx <- t(mtx)
data@x_test <- as.data.frame(mtx)
colnames(data@x_test) <- gsub("-", ".", colnames(data@x_test))
gc()

predictions <- predict(sc_ElasticNet_FS_results, newx = as.matrix(data@x_test), s = "lambda.1se", type = "class")
predictions <- as.factor(ifelse(predictions[, 1] == 0, "ZZZ", snakemake@config[["label_of_interest"]]))
cat("Confusion matrix of ElasticNet predictions on Test1:", "\n")
cm <- confusionMatrix(predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_ElasticNet_Test1 <- cm$overall[["Accuracy"]]
F1_ElasticNet_Test1 <- cm$byClass[["F1"]]

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

genes_with_dot <- setdiff(gsub(".", "-", input_genes, fixed = TRUE), rownames(seu_test2@assays$RNA@data))
mtx_genes <- intersect(gsub(".", "-", input_genes, fixed = TRUE), rownames(seu_test2@assays$RNA@data))
mtx_genes <- c(mtx_genes, gsub("-", ".", genes_with_dot))
mtx <- seu_test2@assays$RNA@data[mtx_genes, ]
mtx <- t(mtx)
x_test2 <- as.data.frame(mtx)
colnames(x_test2) <- gsub("-", ".", colnames(x_test2))
y_test2 <- seu_test2$Diagnosis
y_test2 <- ifelse(y_test2 == snakemake@config[["label_of_interest"]], snakemake@config[["label_of_interest"]], "ZZZ")
y_test2 <- factor(y_test2)

predictions <- predict(sc_ElasticNet_FS_results, newx = as.matrix(x_test2), s = "lambda.1se", type = "class")
predictions <- as.factor(ifelse(predictions[, 1] == 0, "ZZZ", snakemake@config[["label_of_interest"]]))
cat("Confusion matrix of ElasticNet predictions on Test2:", "\n")
cm <- confusionMatrix(predictions, y_test2, mode = "everything")  
print(cm)
Accuracy_ElasticNet_Test2 <- cm$overall[["Accuracy"]]
F1_ElasticNet_Test2 <- cm$byClass[["F1"]]


### Train and test RF model with selected features
rf_model <- ranger(y ~ ., data = xy_train[, c(names(selected_features_ElasticNet), "y")], importance = "impurity")

predictions <- predict(rf_model, data = data@x_test, type = "response")$predictions
cat("Confusion matrix of RF model predictions on Test1:", "\n")
cm <- confusionMatrix(predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_RF_Test1 <- cm$overall[["Accuracy"]]
F1_RF_Test1 <- cm$byClass[["F1"]]

predictions <- predict(rf_model, data = x_test2, type = "response")$predictions
cat("Confusion matrix of RF model predictions on Test2:", "\n")
cm <- confusionMatrix(predictions, y_test2, mode = "everything")  
print(cm)
Accuracy_RF_Test2 <- cm$overall[["Accuracy"]]
F1_RF_Test2 <- cm$byClass[["F1"]]

### Train and test kNN classifier for comparison across feature selection methods
set.seed(123)
train_control <- trainControl(method = "cv", number = 10)
cl <- makePSOCKcluster(snakemake@threads)
registerDoParallel(cl)
start_time <- Sys.time()
knn_model <- train(y ~ ., data = xy_train[, c(names(selected_features_ElasticNet), "y")], 
                   method = "knn", trControl = train_control,
                   tuneLength = 10)
end_time <- Sys.time()
stopCluster(cl)
print("kNN training finished after:")
print(end_time - start_time)
print(knn_model)

knn_predictions <- predict(knn_model, newdata = data@x_test)
cat("Confusion matrix of ElasticNet features kNN model on Test1:", "\n")
cm <- confusionMatrix(knn_predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_kNN_Test1 <- cm$overall[["Accuracy"]]
F1_kNN_Test1 <- cm$byClass[["F1"]]

knn_predictions <- predict(knn_model, newdata = x_test2)
cat("Confusion matrix of ElasticNet features kNN model on Test2:", "\n")
cm <- confusionMatrix(knn_predictions, y_test2, mode = "everything")  
print(cm)
Accuracy_kNN_Test2 <- cm$overall[["Accuracy"]]
F1_kNN_Test2 <- cm$byClass[["F1"]]

### Save results
elasticnet_results <- data.frame(gene = names(selected_features_ElasticNet), coef = selected_features_ElasticNet)
elasticnet_results$Accuracy_ElasticNet_Test1 <- Accuracy_ElasticNet_Test1
elasticnet_results$Accuracy_ElasticNet_Test2 <- Accuracy_ElasticNet_Test2
elasticnet_results$Accuracy_kNN_Test1 <- Accuracy_kNN_Test1
elasticnet_results$Accuracy_kNN_Test2 <- Accuracy_kNN_Test2
elasticnet_results$F1_ElasticNet_Test1 <- F1_ElasticNet_Test1
elasticnet_results$F1_ElasticNet_Test2 <- F1_ElasticNet_Test2
elasticnet_results$F1_kNN_Test1 <- F1_kNN_Test1
elasticnet_results$F1_kNN_Test2 <- F1_kNN_Test2
elasticnet_results$Accuracy_RF_Test1 <- Accuracy_RF_Test1
elasticnet_results$Accuracy_RF_Test2 <- Accuracy_RF_Test2
elasticnet_results$F1_RF_Test1 <- F1_RF_Test1
elasticnet_results$F1_RF_Test2 <- F1_RF_Test2

write.csv(elasticnet_results, file = snakemake@output[[1]], row.names = FALSE)