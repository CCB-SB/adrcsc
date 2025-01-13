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

### scGeneRanger with hyperparameter and input feature set tuning
sc_RF_RFE_with_tuning_grid <- function(x, y, x_v2, y_v2, hyper_grid, input_feature_list, nfolds = 5, step = 0.2) {
  
  set.seed(42)
  results_list <- list()
  data <- data.frame(cbind(x, y = y))
  runtimes <- c()
  acc_on_V2 <- c()
  set_size <- c()
  colnames(x_v2) <- gsub("-", ".", colnames(x_v2))
  colnames(data) <- gsub("-", ".", colnames(data))

  print("Dimensions of trainV1:")
  print(dim(x))
  print("Dimensions of V2 (validation):")
  print(dim(x_v2))
  
  for (i in 1:nrow(hyper_grid)) {
    # Extract hyperparameter configuration
    ntree <- hyper_grid[i, "ntree"]
    maxdepth <- hyper_grid[i, "maxdepth"]
    feature_set <- input_feature_list[[hyper_grid[i, "feature_set"]]]
    feature_set <- gsub("-", ".", feature_set)
    
    cat("Running scGeneRanger with ntree = ", ntree, ", maxdepth = ", maxdepth, ", input feature set = ", hyper_grid[i, "feature_set"], "\n")
  
    print(setdiff(c(feature_set, "y"), colnames(data)))
    feature_set <- intersect(feature_set, colnames(data))

    # Define custom rfe functions for ranger
    if (maxdepth > 0) {
      rangerFuncs <-  list(
        summary = defaultSummary,
        fit = function(x, y, first, last, ...) {
          loadNamespace("ranger")
          data <- data.frame(cbind(x, y = as.factor(y)))
          model <- ranger::ranger(y ~ .,
                                  data = data,
                                  importance = "impurity",
                                  num.trees = ntree,
                                  max.depth = maxdepth)
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
                                  num.trees = ntree)
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

    sizes <- convert_stepsize(length(feature_set), step, min_size = 3)

    # Perform feature selection
    rfe_control <- rfeControl(functions = rangerFuncs, method = "cv", number = nfolds, 
                              returnResamp = "all", saveDetails = TRUE,
                              verbose = TRUE)

    cl <- makePSOCKcluster(snakemake@threads)
    registerDoParallel(cl)

    start_time <- Sys.time()
    results <- rfe(form = as.formula(y ~ .), data = data[, c(feature_set, "y")], sizes = sizes, rfeControl = rfe_control)
    end_time <- Sys.time()

    if (results$optsize > 500) {
      size <- getMinimalSubsetWithinTolerance(results, tolerance = 0.01)
      selected_features <- getAverageFeatureImps(results$variables, size)
    } else if (results$optsize > 150) {
      size <- getMinimalSubsetWithinTolerance(results, tolerance = 0.005)
      selected_features <- getAverageFeatureImps(results$variables, size)
    } else {
      selected_features <- getAverageFeatureImps(results$variables, results$optsize)
    }

    # Train final model on full TrainV1 with selected features
    results$fit <- ranger(y ~ ., data = data[, c(selected_features$var, "y")], importance = "impurity")

    stopCluster(cl)

    # Assess performance on validation set V2
    predictions <- predict(results$fit, data = x_v2[, selected_features$var], type = "response")$predictions
    cat("Confusion matrix of final model on V2:", "\n")
    cm <- confusionMatrix(predictions, y_v2, mode = "everything")  
    print(cm)
    accuracy <- cm$overall[["Accuracy"]]
    
    # Store results 
    results$accuracy_on_V2 <- accuracy
    results$hyperparameter_config <- c(ntree, maxdepth, hyper_grid[i, "feature_set"])
    results$selected_genes <- selected_features$var
    results$importance_scores <- selected_features$Importance
    results$runtime <- end_time - start_time

    results_list[[i]] <- results
    runtimes <- c(runtimes, end_time - start_time)
    acc_on_V2 <- c(acc_on_V2, accuracy)
    set_size <- c(set_size, length(results$selected_genes))
  }
  hyper_grid$accuracy_on_V2 <- acc_on_V2
  hyper_grid$runtime <- runtimes
  hyper_grid$set_size <- set_size
  write.csv(hyper_grid, file = snakemake@output[[2]])
  return(results_list)
}



### Load and subset data
load(snakemake@input[[1]])
seu <- readRDS(snakemake@input[[2]])
seu_trainV1 <- subset(seu, cells = rownames(data@x_train))
DefaultAssay(seu_trainV1) <- "RNA"
seu_trainV1 <- NormalizeData(seu_trainV1)
seu_V2 <- subset(seu, cells = rownames(data@x_v2))

### Define input feature sets
seu_trainV1 <- FindVariableFeatures(seu_trainV1, selection.method = "vst", nfeatures = 2000)
var_features_top2000 <- VariableFeatures(seu_trainV1)
seu_trainV1 <- FindVariableFeatures(seu_trainV1, selection.method = "vst", nfeatures = 4000)
var_features_top4000 <- VariableFeatures(seu_trainV1)
print(length(var_features_top4000))
gene_names <- rownames(seu_trainV1)
l <- lapply(1:3, function(x) sample(gene_names, 2000, replace = FALSE))
input_features <- list()
input_features[[1]] <- l[[1]]
input_features[[2]] <- l[[2]]
input_features[[3]] <- l[[3]]
input_features[[4]] <- var_features_top2000
input_features[[5]] <- var_features_top4000

### Get expression values for all input feature sets
all_input_features <- unique(c(input_features[[1]], input_features[[2]], input_features[[3]],
                               input_features[[4]], input_features[[5]]))
mtx <- seu_trainV1@assays$RNA@data[all_input_features, ]
mtx <- t(mtx)
data@x_train <- as.data.frame(mtx)
mtx <- seu_V2@assays$RNA@data[all_input_features, ]
mtx <- t(mtx)
data@x_v2 <- as.data.frame(mtx)
gc()

### Define hyperparameter grid
hyperparameter_grid <- expand.grid(
  ntree = c(100, 500, 1000),      # Number of trees in the forest
  maxdepth = c(10, 20, 0),        # Maximum depth of trees
  feature_set = c("random2000_1", "random2000_2", "random2000_3", "top2000HVGs", "top4000HVGs")  # Input feature set names
)

### Call scGeneRanger with hyperparameter grid
scGR_results <- sc_RF_RFE_with_tuning_grid(x=data@x_train, y=data@y_train, x_v2=data@x_v2, y_v2=data@y_v2, hyper_grid=hyperparameter_grid, input_feature_list=input_features)
save(scGR_results, file = snakemake@output[[1]])
