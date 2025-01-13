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
  threshold <- rfe_results$results$WeightedF1[rfe_results$results$Variables == rfe_results$optsize]  - tolerance
  results_within_tolerance <- rfe_results$results[rfe_results$results$WeightedF1 >= threshold, ]
  min_subset <- results_within_tolerance$Variables[which.min(results_within_tolerance$Variables)]
  min_subset
}

### F1-score summary function with weighting of validation fold sizes
weighted_f1 <- function(data, lev = NULL, model = NULL, fold_sizes = NULL, nfolds) {
  # Calculate precision, recall and F1-score
  positive <- lev[1]
  precision <- posPredValue(data$pred, data$obs, positive = positive)
  recall <- sensitivity(data$pred, data$obs, positive = positive)
  f1_val <- (2 * precision * recall) / (precision + recall)
  
  # Weight the F1-score by fold sizes (and multiply by nfolds to offset averaging across folds by rfe())
  fold_weight <- length(data$obs) / sum(fold_sizes)
  weighted_f1 <- f1_val * fold_weight * nfolds
  
  names(weighted_f1) <- "WeightedF1"
  return(weighted_f1)
}

### scGeneRanger
sc_RF_RFE <- function(x, y, num.trees, max.depth, folds, step = 0.2) {
  
  set.seed(42)
  nfolds <- length(folds)
  train_indices <- lapply(folds, function(x) x$train)
  valid_indices <- lapply(folds, function(x) x$valid)
  fold_sizes <- sapply(folds, function(x) length(x$valid))
  data <- data.frame(cbind(x, y = y))
  colnames(data) <- gsub("-", ".", colnames(data))

  print("Dimensions of complete training set:")
  print(dim(x))
  
  cat("Running scGeneRanger with ntree = ", num.trees, "maxdepth = ", max.depth, "\n")

  # Define custom rfe functions for ranger
  if (max.depth == 0) {
    rangerFuncs <-  list(
      summary = function(data, lev = NULL, model = NULL) {
        weighted_f1(data, lev, model, fold_sizes, nfolds)
      },
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
      summary = function(data, lev = NULL, model = NULL) {
        weighted_f1(data, lev, model, fold_sizes, nfolds)
      },
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
  seeds <- lapply(1:nfolds, function(x) sample.int(1000, length(sizes)))
  seeds[[nfolds+1]] <- sample.int(1000, 1)

  # Perform feature selection
  rfe_control <- rfeControl(functions = rangerFuncs, method = "cv", number = nfolds, 
                            returnResamp = "all", saveDetails = TRUE, allowParallel = TRUE,
                            verbose = TRUE, seeds = seeds, index = train_indices, indexOut = valid_indices)

  cl <- makePSOCKcluster(snakemake@threads)
  clusterExport(cl, "weighted_f1")
  registerDoParallel(cl)

  start_time <- Sys.time()
  results <- rfe(form = as.formula(y ~ .), data = data, sizes = sizes, rfeControl = rfe_control, metric = "WeightedF1")
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

load(snakemake@input[["fold_indices"]])
train_rownames <- unique(unlist(lapply(fold_indices, function(x) c(x$train, x$valid))))
seu_train <- subset(seu_train, cells = train_rownames)
seu_train

### Extract best hyperparam config
print(length(snakemake@config[["celltypes"]]))
hyperparam_results <- data.frame()
for (i in 1:length(snakemake@config[["celltypes"]])) {
  hyperparam_df <- rbind(read.csv(snakemake@input[["hyperparams1"]][[i]]), read.csv(snakemake@input[["hyperparams2"]][[i]]))
  celltype <- str_extract(snakemake@input[["hyperparams1"]][[i]], "(?<=results/)[^/]+")
  hyperparam_df$celltype <- celltype
  hyperparam_results <- rbind(hyperparam_results, hyperparam_df)
}

print(head(hyperparam_results))
print(tail(hyperparam_results))

best_config <- hyperparam_results %>% 
                group_by(ntree, maxdepth, feature_set) %>%
                summarize(mean_F1_on_V2 = mean(F1_on_V2)) %>%
                arrange(desc(mean_F1_on_V2)) %>%
                slice(1)

#write.csv(best_config, file = "hypparam_best_config.csv")

ntree <- 500 #best_config$ntree[1] 
maxdepth <- 20 #best_config$maxdepth[1]
feature_set <- "LowExpressionFilter" #best_config$feature_set[1] #best performing configuration in the hyperparameter optimization included a randomly chosen input feature set, 
                                                                 #so the best configuration not including a randomly chosen input feature set was manually chosen here 

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
  print("oh oh. no good.")
  print(feature_set)
  stop("Random set performed best :(")
}

### Get expression values and labels
mtx <- seu_train@assays$RNA@data[input_features, ]
mtx <- t(mtx)
data@x_train <- as.data.frame(mtx)

y <- seu_train$Diagnosis
y <- ifelse(y == snakemake@config[["label_of_interest"]], snakemake@config[["label_of_interest"]], "ZZZ")
data@y_train <- factor(y)

gc()

### Feature selection
scGR_results <- sc_RF_RFE(x=data@x_train, y=data@y_train, num.trees = ntree, max.depth = maxdepth, folds = fold_indices)
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
set.seed(42)
xy_train <- cbind(data@x_train, y = data@y_train)
balanced_train_data <- downSample(xy_train, data@y_train, yname = "y")
dim(balanced_train_data)
colnames(balanced_train_data) <- gsub("-", ".", colnames(balanced_train_data))
scGR_results$fit <- ranger(y ~ ., data = balanced_train_data[, c(selected_features$var, "y")], importance = "impurity")

### Assess performance on test sets Test1 and Test2
seu_test1 <- subset(seu, cells = rownames(data@x_test))
seu_test2 <- readRDS(file = snakemake@input[[3]])
DefaultAssay(seu_test1) <- "RNA"
DefaultAssay(seu_test2) <- "RNA"
seu_test1 <- NormalizeData(seu_test1)
seu_test2 <- NormalizeData(seu_test2)

# Test1
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

# Test2
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
knn_model <- train(y ~ ., data = balanced_train_data[, c(selected_features$var, "y")], 
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
