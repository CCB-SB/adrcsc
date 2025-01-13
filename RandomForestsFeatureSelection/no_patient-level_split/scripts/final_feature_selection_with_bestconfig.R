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
suppressMessages(library(stringr))

set.seed(123)

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile, type = "output")
sink(logfile, type = "message")

### Custom data class
DataClass <- setClass("DataClass", slots = c(x_train = "data.frame", y_train = "factor", x_v2 = "data.frame", y_v2 = "factor", x_test = "data.frame", y_test = "factor"))

### Helper function for rfe
convert_stepsize <- function(n, step, min_size = 5) {
  sizes <- c()
  current_index <- n
  while (current_index >= min_size) {
    sizes <- c(sizes, current_index)
    current_index <- current_index * (1 - step)
  }
  return(round(sizes))
}

### Helper function that extracts the average feature importances from the rfe object
getAverageFeatureImps <- function(vardf, size) {
  finalImp <- ddply(vardf[, c("Overall", "var")],
                    .(var),
                    function(x) mean(x$Overall, na.rm = TRUE))
  names(finalImp)[2] <- "Importance"
  finalImp <- finalImp[order(finalImp$Importance, decreasing = TRUE),]
  finalImp <- finalImp[1:size,]
  finalImp
}

### Helper function
getMinimalSubsetWithinTolerance <- function(rfe_results, tolerance = 0.005) {
  accuracy_threshold <- rfe_results$results$Accuracy[rfe_results$results$Variables == rfe_results$optsize]  - tolerance
  results_within_tolerance <- rfe_results$results[rfe_results$results$Accuracy >= accuracy_threshold, ]
  min_subset <- results_within_tolerance$Variables[which.min(results_within_tolerance$Variables)]
  min_subset
}

### scGeneRanger
sc_RF_RFE <- function(x, y, num.trees, max.depth, nfolds = 5, step = 0.2) {
  
  set.seed(42)
  data <- data.frame(cbind(x, y = y))
  colnames(data) <- gsub("-", ".", colnames(data))

  print("Dimensions of complete training set:")
  print(dim(x))
  
  cat("Running scGeneRanger with ntree = ", num.trees, "maxdepth = ", max.depth, "\n")

  # Define custom rfe functions for ranger
  if (max.depth == 0) {
    rangerFuncs <-  list(
      summary = defaultSummary,
      fit = function(x, y, first, last, ...) {
        loadNamespace("ranger")
        data <- data.frame(cbind(x, y = as.factor(y)))
        model <- ranger::ranger(y ~ .,
                                data = data,
                                importance = "impurity",
                                num.trees = num.trees)
        return(model)
      },
      pred = function(object, x)  {
        predict(object$forest, data = x, type = "response")$predictions
      },
      rank = function(object, x, y) {
        vimp <- data.frame(Overall = object$variable.importance)
        vimp <- vimp[order(vimp$Overall, decreasing = TRUE),, drop = FALSE]
        if (ncol(x) == 1) {
          vimp$var <- colnames(x)
        } else {
          vimp$var <- rownames(vimp)
        }
        vimp
      },
      selectSize = pickSizeBest,
      selectVar = pickVars
    )
  } else {
    rangerFuncs <-  list(
      summary = defaultSummary,
      fit = function(x, y, first, last, ...) {
        loadNamespace("ranger")
        data <- data.frame(cbind(x, y = as.factor(y)))
        model <- ranger::ranger(y ~ .,
                                data = data,
                                importance = "impurity",
                                num.trees = num.trees,
                                max.depth = max.depth)
        return(model)
      },
      pred = function(object, x)  {
        predict(object$forest, data = x, type = "response")$predictions
      },
      rank = function(object, x, y) {
        vimp <- data.frame(Overall = object$variable.importance)
        vimp <- vimp[order(vimp$Overall, decreasing = TRUE),, drop = FALSE]
        if (ncol(x) == 1) {
          vimp$var <- colnames(x)
        } else {
          vimp$var <- rownames(vimp)
        }
        vimp
      },
      selectSize = pickSizeBest,
      selectVar = pickVars
    ) 
  }

  sizes <- convert_stepsize(dim(x)[2], step, min_size = 3)

  # Perform feature selection
  rfe_control <- rfeControl(functions = rangerFuncs, method = "cv", number = nfolds, 
                            returnResamp = "all", saveDetails = TRUE,
                            verbose = TRUE)

  cl <- makePSOCKcluster(snakemake@threads)
  registerDoParallel(cl)

  start_time <- Sys.time()
  results <- rfe(form = as.formula(y ~ .), data = data, sizes = sizes, rfeControl = rfe_control)
  end_time <- Sys.time()

  print("scGeneRanger finished after:")
  print(end_time - start_time)

  stopCluster(cl)
  return(results)
}



### Load and subset data
load(snakemake@input[[1]])
seu <- readRDS(snakemake@input[[2]])
seu_train <- subset(seu, cells = c(rownames(data@x_train), rownames(data@x_v2)))
DefaultAssay(seu_train) <- "RNA"
seu_train <- NormalizeData(seu_train)

### Extract best hyperparam config
hyperparam_results <- read.csv(file = snakemake@input[["hyperparams"]])

ntree <- hyperparam_results$ntree[1]
maxdepth <- hyperparam_results$maxdepth[1]
feature_set <- hyperparam_results$feature_set[1]

cat("Hyperparameter config with best performance on validation set V2:", "\n")
cat("ntree: ", ntree, "maxdepth: ", maxdepth, "feature_set: ", feature_set, "\n")

if (feature_set == "top2000HVGs") {
  seu_train <- FindVariableFeatures(seu_train, selection.method = "vst", nfeatures = 2000)
  input_features <- VariableFeatures(seu_train)
} else if (feature_set == "top4000HVGs") {
  seu_train <- FindVariableFeatures(seu_train, selection.method = "vst", nfeatures = 4000)
  input_features <- VariableFeatures(seu_train)
} else if (feature_set == "top6000HVGs") {
  seu_train <- FindVariableFeatures(seu_train, selection.method = "vst", nfeatures = 6000)
  input_features <- VariableFeatures(seu_train)
} else if (feature_set == "top8000HVGs") {
  seu_train <- FindVariableFeatures(seu_train, selection.method = "vst", nfeatures = 8000)
  input_features <- VariableFeatures(seu_train)
} else if (feature_set == "LowExpressionFilter"){
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
  input_features <- names(genes_to_keep[genes_to_keep])
} else if (feature_set == "DEGs") {
  ### Load DEGs table and get sig. DEGs for current celltype
  DEGs <- read.csv(file = snakemake@input[[1]])
  celltype <- snakemake@config[["celltypes"]][[snakemake@params[[1]]]][["name"]]
  DEGs <- DEGs %>% filter(cluster_id == celltype, category != "Not deregulated")
  input_features <- DEGs$gene
} else {
  print("A random input feature set resulted in the highest performance: ")
  print(feature_set)
}

### Get expression values and labels
mtx <- seu_train@assays$RNA@data[input_features, ]
mtx <- t(mtx)
data@x_train <- as.data.frame(mtx)

common_levels <- union(levels(data@y_train), levels(data@y_v2))
y_trainV1V2 <- factor(c(as.character(data@y_train), as.character(data@y_v2)), levels = common_levels)
data@y_train <- y_trainV1V2

gc()

### Feature selection
scGR_results <- sc_RF_RFE(x=data@x_train, y=data@y_train, num.trees = ntree, max.depth = maxdepth)
if (scGR_results$optsize > 500) {
  size <- getMinimalSubsetWithinTolerance(scGR_results, tolerance = 0.01)
  selected_features <- getAverageFeatureImps(scGR_results$variables, size)
} else if (scGR_results$optsize > 150) {
  size <- getMinimalSubsetWithinTolerance(scGR_results, tolerance = 0.005)
  selected_features <- getAverageFeatureImps(scGR_results$variables, size)
} else {
  selected_features <- getAverageFeatureImps(scGR_results$variables, scGR_results$optsize)
}

# Save backwards feature elimination cv performance plot
subset_perf_plot <- plot(scGR_results, type = c("g", "o"))
ggsave(file.path(snakemake@output[[4]]), plot = as.ggplot(subset_perf_plot), width = 6, height = 4, dpi = 300)

### Train final model on full TrainV1V2 with final feature selection
xy_train <- cbind(data@x_train, y = data@y_train)
colnames(xy_train) <- gsub("-", ".", colnames(xy_train))
scGR_results$fit <- ranger(y ~ ., data = xy_train[, c(selected_features$var, "y")], importance = "impurity")

### Assess performance on test sets Test1 and Test2
seu_test1 <- subset(seu, cells = rownames(data@x_test))
seu_test2 <- readRDS(file = snakemake@input[[3]])
DefaultAssay(seu_test1) <- "RNA"
DefaultAssay(seu_test2) <- "RNA"
seu_test1 <- NormalizeData(seu_test1)
seu_test2 <- NormalizeData(seu_test2)

# Test1 (Test cells from same patients)
genes_with_dot <- setdiff(gsub(".", "-", selected_features$var, fixed = TRUE), rownames(seu_test1@assays$RNA@data))
mtx_genes <- intersect(gsub(".", "-", selected_features$var, fixed = TRUE), rownames(seu_test1@assays$RNA@data))
mtx_genes <- c(mtx_genes, gsub("-", ".", genes_with_dot))
mtx <- seu_test1@assays$RNA@data[mtx_genes, ]
mtx <- t(mtx)
data@x_test <- as.data.frame(mtx)
colnames(data@x_test) <- gsub("-", ".", colnames(data@x_test))
gc()

predictions <- predict(scGR_results$fit, data = data@x_test, type = "response")$predictions
cat("Confusion matrix of final RF model on Test1:", "\n")
cm <- confusionMatrix(predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_RF_Test1 <- cm$overall[["Accuracy"]]
F1_RF_Test1 <- cm$byClass[["F1"]]

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

genes_with_dot <- setdiff(gsub(".", "-", selected_features$var, fixed = TRUE), rownames(seu_test2@assays$RNA@data))
mtx_genes <- intersect(gsub(".", "-", selected_features$var, fixed = TRUE), rownames(seu_test2@assays$RNA@data))
mtx_genes <- c(mtx_genes, gsub("-", ".", genes_with_dot))
mtx <- seu_test2@assays$RNA@data[mtx_genes, ]
mtx <- t(mtx)
x_test2 <- as.data.frame(mtx)
colnames(x_test2) <- gsub("-", ".", colnames(x_test2))
y_test2 <- seu_test2$Diagnosis
y_test2 <- ifelse(y_test2 == snakemake@config[["label_of_interest"]], snakemake@config[["label_of_interest"]], "ZZZ")
y_test2 <- factor(y_test2)

predictions <- predict(scGR_results$fit, data = x_test2, type = "response")$predictions
cat("Confusion matrix of final RF model on Test2:", "\n")
cm <- confusionMatrix(predictions, y_test2, mode = "everything")  
print(cm)
Accuracy_RF_Test2 <- cm$overall[["Accuracy"]]
F1_RF_Test2 <- cm$byClass[["F1"]]

### Train kNN classifier for comparison across feature selection methods
set.seed(42)
train_control <- trainControl(method = "cv", number = 10)
cl <- makePSOCKcluster(snakemake@threads)
registerDoParallel(cl)
start_time <- Sys.time()
knn_model <- train(y ~ ., data = xy_train[, c(selected_features$var, "y")], 
                   method = "knn", trControl = train_control,
                   tuneLength = 10)
end_time <- Sys.time()
stopCluster(cl)
print("kNN training finished after:")
print(end_time - start_time)
print(knn_model)
png(file = snakemake@output[[3]])
plot(knn_model)
dev.off()

knn_predictions <- predict(knn_model, newdata = data@x_test)
cat("Confusion matrix of kNN model on Test1:", "\n")
cm <- confusionMatrix(knn_predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_kNN_Test1 <- cm$overall[["Accuracy"]]
F1_kNN_Test1 <- cm$byClass[["F1"]]

knn_predictions <- predict(knn_model, newdata = x_test2)
cat("Confusion matrix of kNN model on Test2:", "\n")
cm <- confusionMatrix(knn_predictions, y_test2, mode = "everything")  
print(cm)
Accuracy_kNN_Test2 <- cm$overall[["Accuracy"]]
F1_kNN_Test2 <- cm$byClass[["F1"]]

### Save results
results_df <- as.data.frame(cbind(selected_features$var, selected_features$Importance))
colnames(results_df) <- c("gene", "imp")
results_df$Accuracy_RF_Test1 <- Accuracy_RF_Test1
results_df$Accuracy_RF_Test2 <- Accuracy_RF_Test2
results_df$Accuracy_kNN_Test1 <- Accuracy_kNN_Test1
results_df$Accuracy_kNN_Test2 <- Accuracy_kNN_Test2
results_df$F1_RF_Test1 <- F1_RF_Test1
results_df$F1_RF_Test2 <- F1_RF_Test2
results_df$F1_kNN_Test1 <- F1_kNN_Test1
results_df$F1_kNN_Test2 <- F1_kNN_Test2

write.csv(results_df, file = snakemake@output[[1]], row.names = FALSE)
save(scGR_results, file = snakemake@output[[2]])
