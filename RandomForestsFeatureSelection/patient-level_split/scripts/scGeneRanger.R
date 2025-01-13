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

### Helper function that extracts the gene expression data, response variable and sample/donor information from a seurat object
ExtractData <- function(seu, genes, response_var, sample_col) {
  data <- seu@assays$RNA@data[genes, ]
  data <- t(data) 
  gc()
  y <- seu[[response_var]]
  sample <- seu[[sample_col]]
  data <- cbind(data, y)
  data <- cbind(data, sample)
  data[, response_var] <- factor(data[, response_var])
  data
}

### Downsampling function for one-vs-all/binary classification
DownsampleIndicesForOneVsAll <- function(indices, data, label_of_interest, n) {
  label_counts <- table(data[indices, ])
  labels <- names(sort(label_counts, decreasing = FALSE))
  downsampled_indices <- integer(0)
  selection_size <- round(n / (length(labels)-1))
  i <- 1
  for (label in labels) {
    label_indices <- indices[data[indices, ] == label]
    if (label == label_of_interest) {
      selected_indices <- sample(label_indices, n)
    } else {
      if (length(label_indices) >= selection_size) {
        selected_indices <- sample(label_indices, selection_size)
      } else {
        selected_indices <- label_indices
        selection_size <- selection_size + (selection_size - length(label_indices)) / (length(labels)-1 - i)
        i <- i+1
      }
    }
    downsampled_indices <- c(downsampled_indices, selected_indices)
  }
  #print(table(data[downsampled_indices, ]))
  downsampled_indices
}

### Helper function that sets binary labels for one-vs-all classification
SetBinaryLabels <- function(data, class_of_interest, label_of_interest) {
  data$binary_label <- "ZZZ"
  data[data[, class_of_interest] == label_of_interest, "binary_label"] <- label_of_interest
  data$binary_label <- factor(data$binary_label)
  data
}

### Create group cv folds based on samples/donors and balance classes in each fold via downsampling
CreateGroupFoldIndices <- function(seu, sample_col, class_of_interest, label_of_interest, k) {
  # Get the unique sample/donor IDs for both classes
  neg_class_samples <- unique(seu[[sample_col]][, 1][seu[[class_of_interest]][, 1] != label_of_interest])
  pos_class_samples <- unique(seu[[sample_col]][, 1][seu[[class_of_interest]][, 1] == label_of_interest])
  
  # Calculate the fold sizes and the number of remaining samples
  fold_size_neg <- floor(length(neg_class_samples) / k)
  fold_size_pos <- floor(length(pos_class_samples) / k)

  remainder_neg <- length(neg_class_samples) %% k
  remainder_pos <- length(pos_class_samples) %% k
  
  # Create a list to store the fold indices
  fold_indices <- list()
  
  # Create k folds
  for (i in 1:k) {
    # Get the indices for the current fold
    neg_idx <- ((i-1) * fold_size_neg + 1):min(i * fold_size_neg, length(neg_class_samples))
    pos_idx <- ((i-1) * fold_size_pos + 1):min(i * fold_size_pos, length(pos_class_samples))
    
    # Include remainders by adding one more sample from the remainder pool to the first few folds
    if (i <= remainder_neg) {
      neg_idx <- c(neg_idx, fold_size_neg * k + i)
    }
    if (i <= remainder_pos) {
      pos_idx <- c(pos_idx, fold_size_pos * k + i)
    }
    
    validation_samples <- c(neg_class_samples[neg_idx], pos_class_samples[pos_idx])
    training_samples <- setdiff(c(neg_class_samples, pos_class_samples), validation_samples)
    
    # Get indices for training and validation
    training_indices <- rownames(seu[[]][seu[[]][[sample_col]] %in% training_samples, ])
    validation_indices <- rownames(seu[[]][seu[[]][[sample_col]] %in% validation_samples, ])
    
    # Downsample the cells to balance the classes in each fold
    n <- min(sum(seu[[class_of_interest]][training_indices, ] == label_of_interest), sum(seu[[class_of_interest]][training_indices, ] != label_of_interest))
    training_indices <- DownsampleIndicesForOneVsAll(training_indices, seu[[class_of_interest]], label_of_interest, n)

    n <- min(sum(seu[[class_of_interest]][validation_indices, ] == label_of_interest), sum(seu[[class_of_interest]][validation_indices, ] != label_of_interest))
    validation_indices <- DownsampleIndicesForOneVsAll(validation_indices, seu[[class_of_interest]], label_of_interest, n)

    # Store the training and validation indices in the fold_indices list
    fold_indices[[i]] <- list("train" = training_indices, "valid" = validation_indices)
    
    # Print fold information
    print(paste("### Fold", i, "###"))

    print("Training set sample/donor IDs:")
    cat(training_samples, sep = ", ")
    cat("\n")
    print(paste0("Number of control samples/donors in training set:"))
    print(sum(training_samples %in% neg_class_samples))
    print(paste0("Number of ", group2, " samples/donors in training set:"))
    print(sum(training_samples %in% pos_class_samples))
    
    print("Size of training set:")
    print(length(training_indices))
    print(paste0("Number of control cells in training set:"))
    print(sum(seu[[]][training_indices, class_of_interest] != label_of_interest))
    print(paste0("Number of ", group2, " cells in training set:"))
    print(sum(seu[[]][training_indices, class_of_interest] == label_of_interest))
    
    print("Validation set sample/donor IDs:")
    cat(paste(validation_samples), sep = ", ")
    cat("\n")
    print(paste0("Number of control samples/donors in validation set:"))
    print(sum(validation_samples %in% neg_class_samples))
    print(paste0("Number of ", group2, " samples/donors in validation set:"))
    print(sum(validation_samples %in% pos_class_samples))
    
    print("Size of validation set:")
    print(length(validation_indices))
    print(paste0("Number of control cells in validation set:"))
    print(sum(seu[[]][validation_indices, class_of_interest] != label_of_interest))
    print(paste0("Number of ", group2, " cells in validation set:"))
    print(sum(seu[[]][validation_indices, class_of_interest] == label_of_interest))
  }
  
  return(fold_indices)
}

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
  threshold <- rfe_results$results$WeightedF1[rfe_results$results$Variables == rfe_results$optsize] - tolerance
  results_within_tolerance <- rfe_results$results[rfe_results$results$WeightedF1 >= threshold, ]
  min_subset <- results_within_tolerance$Variables[which.min(results_within_tolerance$Variables)]
  min_subset
}

### F1-score summary function for rfe()
# f1 <- function(data, lev = NULL, model = NULL) {
#   positive <- lev[1]
#   precision <- posPredValue(data$pred, data$obs, positive = positive)
#   recall <- sensitivity(data$pred, data$obs, positive = positive)
#   f1_val <- (2*precision*recall) / (precision + recall)
#   names(f1_val) <- c("F1")
#   f1_val
# }

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

### scGeneRanger with hyperparameter and input feature set tuning
sc_RF_RFE_with_tuning_grid <- function(x, y, x_v2, y_v2, hyper_grid, input_feature_list, folds, step = 0.2) {
  
  set.seed(42)
  results_list <- list()
  nfolds <- length(folds)
  train_indices <- lapply(folds, function(x) x$train)
  valid_indices <- lapply(folds, function(x) x$valid)
  fold_sizes <- sapply(folds, function(x) length(x$valid))
  data <- data.frame(cbind(x, y = y))
  runtimes <- c()
  F1_on_V2 <- c()
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

    # Define custom rfe functions for scGeneRanger
    if (maxdepth > 0) {
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
        summary = function(data, lev = NULL, model = NULL) {
          weighted_f1(data, lev, model, fold_sizes, nfolds)
        },
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
    results <- rfe(form = as.formula(y ~ .), data = data[, c(feature_set, "y")], sizes = sizes, rfeControl = rfe_control, metric = "WeightedF1")
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

    cat("scGeneRanger finished after", as.numeric(difftime(end_time, start_time, units = "mins")), "mins", "\n")

    # Train final model on full TrainV1 with selected features
    balanced_train_data <- downSample(data, y, yname = "y")
    dim(balanced_train_data)
    results$fit <- ranger(y ~ ., data = balanced_train_data[, c(selected_features$var, "y")], importance = "impurity")

    stopCluster(cl)

    # Assess performance on validation set V2
    predictions <- predict(results$fit, data = x_v2[, selected_features$var], type = "response")$predictions
    cat("Confusion matrix of final model on V2:", "\n")
    cm <- confusionMatrix(predictions, y_v2, mode = "everything")  
    print(cm)
    F1 <- cm$byClass[["F1"]]
    
    # Store results 
    results$F1_on_V2 <- F1
    results$hyperparameter_config <- c(ntree, maxdepth, hyper_grid[i, "feature_set"])
    results$selected_genes <- selected_features$var
    results$importance_scores <- selected_features$Importance
    results$runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))

    results_list[[i]] <- results
    runtimes <- c(runtimes, end_time - start_time)
    F1_on_V2 <- c(F1_on_V2, F1)
    set_size <- c(set_size, length(results$selected_genes))
  }
  hyper_grid$F1_on_V2 <- F1_on_V2
  hyper_grid$runtime <- runtimes
  hyper_grid$set_size <- set_size
  write.csv(hyper_grid, file = snakemake@output[[2]])
  return(results_list)
}

meta_data_column <- snakemake@config[["class_of_interest"]]
donor_slot <- snakemake@config[["donor_slot"]]
group2 <- snakemake@config[["label_of_interest"]]

### Load and subset data
load(snakemake@input[[1]])
seu <- readRDS(snakemake@input[[2]])
seu_trainV1 <- subset(seu, cells = rownames(data@x_train))
DefaultAssay(seu_trainV1) <- "RNA"
seu_trainV1 <- NormalizeData(seu_trainV1)
seu_V2 <- subset(seu, cells = rownames(data@x_v2))

fold_indices <- CreateGroupFoldIndices(seu_trainV1, sample_col = donor_slot, class_of_interest = meta_data_column, 
                                       label_of_interest = group2, k = 5)
train_rownames <- unique(unlist(lapply(fold_indices, function(x) c(x$train, x$valid))))
seu_subset <- subset(seu_trainV1, cells = train_rownames)
seu_subset
print(sum(seu_subset[[meta_data_column]][, 1] == group2))
print(sum(seu_subset[[meta_data_column]][, 1] != group2))

### Define input feature sets
seu_subset <- FindVariableFeatures(seu_subset, selection.method = "vst", nfeatures = 2000)
var_features_top2000 <- VariableFeatures(seu_subset)
seu_subset <- FindVariableFeatures(seu_subset, selection.method = "vst", nfeatures = 4000)
var_features_top4000 <- VariableFeatures(seu_subset)
print(length(var_features_top4000))
gene_names <- rownames(seu_subset)
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
mtx <- seu_subset@assays$RNA@data[all_input_features, ]
mtx <- t(mtx)
data@x_train <- as.data.frame(mtx)
y <- seu_subset$Diagnosis
y <- ifelse(y == snakemake@config[["label_of_interest"]], snakemake@config[["label_of_interest"]], "ZZZ")
y <- factor(y)
data@y_train <- y
mtx <- seu_V2@assays$RNA@data[all_input_features, ]
mtx <- t(mtx)
data@x_v2 <- as.data.frame(mtx)
gc()

print(setdiff(gsub("-", ".", input_features[[5]]), gsub("-", ".", colnames(data@x_train))))

### Define hyperparameter grid
hyperparameter_grid <- expand.grid(
  ntree = c(100, 500, 1000),      # Number of trees in the forest
  maxdepth = c(10, 20, 0),        # Maximum depth of trees
  feature_set = c("random2000_1", "random2000_2", "random2000_3", "top2000HVGs", "top4000HVGs")  # Input feature set names
)

### Call scGeneRanger with hyperparameter grid
scGR_results <- sc_RF_RFE_with_tuning_grid(x=data@x_train, y=data@y_train, x_v2=data@x_v2, y_v2=data@y_v2, hyper_grid=hyperparameter_grid, input_feature_list=input_features, folds = fold_indices)
save(scGR_results, file = snakemake@output[[1]])
save(fold_indices, file = snakemake@output[[3]])
