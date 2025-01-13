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
suppressMessages(library(class))

set.seed(123)

### Set up logging
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile, type = "output")
sink(logfile, type = "message")

### Custom data class
DataClass <- setClass("DataClass", slots = c(x_train = "data.frame", y_train = "factor", x_test = "data.frame", y_test = "factor"))

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
  neg_class_samples <- paste(unique(seu[[sample_col]][, 1][seu[[class_of_interest]][, 1] != label_of_interest]))
  pos_class_samples <- paste(unique(seu[[sample_col]][, 1][seu[[class_of_interest]][, 1] == label_of_interest]))
  
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
    cat(paste(training_samples), sep = ", ")
    cat("\n")
    print(paste0("Number of control samples/donors in training set:"))
    print(sum(training_samples %in% neg_class_samples))
    print(paste0("Number of ", positive_class, " samples/donors in training set:"))
    print(sum(training_samples %in% pos_class_samples))
    
    print("Size of training set:")
    print(length(training_indices))
    print(paste0("Number of control cells in training set:"))
    print(sum(seu[[]][training_indices, class_of_interest] != label_of_interest))
    print(paste0("Number of ", positive_class, " cells in training set:"))
    print(sum(seu[[]][training_indices, class_of_interest] == label_of_interest))
    
    print("Validation set sample/donor IDs:")
    cat(paste(validation_samples), sep = ", ")
    cat("\n")
    print(paste0("Number of control samples/donors in validation set:"))
    print(sum(validation_samples %in% neg_class_samples))
    print(paste0("Number of ", positive_class, " samples/donors in validation set:"))
    print(sum(validation_samples %in% pos_class_samples))
    
    print("Size of validation set:")
    print(length(validation_indices))
    print(paste0("Number of control cells in validation set:"))
    print(sum(seu[[]][validation_indices, class_of_interest] != label_of_interest))
    print(paste0("Number of ", positive_class, " cells in validation set:"))
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
sc_RF_RFE <- function(x, y, folds, step = 0.2) {
  
  set.seed(42)
  nfolds <- length(folds)
  train_indices <- lapply(folds, function(x) x$train)
  valid_indices <- lapply(folds, function(x) x$valid)
  fold_sizes <- sapply(folds, function(x) length(x$valid))
  data <- data.frame(cbind(x, y = y))
  colnames(data) <- gsub("-", ".", colnames(data))

  print("Dimensions of complete training set:")
  print(dim(x))
  
  cat("Running scGeneRanger with ntree = 500, maxdepth = NULL, input feature set = top2000HVGs", "\n")

  # Define custom rfe functions for ranger
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
                              num.trees = 500)
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

  sizes <- convert_stepsize(dim(x)[2], step, min_size = 3)
  
  seeds <- lapply(1:nfolds, function(x) sample.int(1000, length(sizes)))
  seeds[[nfolds+1]] <- sample.int(1000, 1)

  # Perform feature selection
  rfe_control <- rfeControl(functions = rangerFuncs, method = "cv", number = nfolds, 
                            returnResamp = "all", saveDetails = TRUE, allowParallel = TRUE,
                            verbose = TRUE, seeds = seeds, index = train_indices, indexOut = valid_indices)

  cl <- makePSOCKcluster(5)
  clusterExport(cl, "weighted_f1")
  registerDoParallel(cl)

  start_time <- Sys.time()
  results <- rfe(form = as.formula(y ~ .), data = data[,], sizes = sizes, rfeControl = rfe_control, metric = "WeightedF1")
  end_time <- Sys.time()

  print("scGeneRanger finished after:")
  print(end_time - start_time)

  stopCluster(cl)
  return(results)
}



### Fetch values from config
class_column <- snakemake@config[["class_column"]]
annotation_column <- snakemake@config[["annotation_column"]]
donor_column <- snakemake@config[["donor_column"]]
positive_class <- snakemake@config[["positive_class"]]
celltype <- snakemake@config[["celltypes"]][[snakemake@params[[1]]]][["name"]]

### Load and subset data
seu_train <- readRDS(file = snakemake@input[[1]])
seu_train
DefaultAssay(seu_train) <- "RNA"
seu_train <- NormalizeData(seu_train)
expr <- FetchData(object = seu_train, vars = annotation_column)
seu_train <- seu_train[, which(x = expr == celltype)]
message(paste(celltype, "Seurat object for training:"))
seu_train
print(table(seu_train[[class_column]]))

seu_test <- readRDS(file = snakemake@input[[2]])
DefaultAssay(seu_test) <- "RNA"
seu_test <- NormalizeData(seu_test)
expr <- FetchData(object = seu_test, vars = annotation_column)
seu_test <- seu_test[, which(x = expr == celltype)]
message(paste(celltype, "Seurat object for independent testing:"))
seu_test
print(table(seu_test[[class_column]]))

### Determine balanced group k-fold CV fold indices based on donors
fold_indices <- CreateGroupFoldIndices(seu_train, sample_col = donor_column, class_of_interest = class_column, 
                                       label_of_interest = positive_class, k = 5)
train_rownames <- unique(unlist(lapply(fold_indices, function(x) c(x$train, x$valid))))
seu_subset <- subset(seu_train, cells = train_rownames)
save(fold_indices, file = snakemake@output[["fold_indices"]])

### Define input feature set (top2000HVGs)
seu_train <- FindVariableFeatures(seu_train, selection.method = "vst", nfeatures = 2000)
input_features <- VariableFeatures(seu_train)

### Get expression matrices and response labels for train and test set
mtx <- seu_subset@assays$RNA@data[input_features, ]
mtx <- t(mtx)
x_train <- as.data.frame(mtx)
y <- seu_subset[[class_column]]
y <- ifelse(y == positive_class, positive_class, "ZZZ")
y_train <- factor(y)

mtx <- seu_test@assays$RNA@data[input_features, ]
mtx <- t(mtx)
x_test <- as.data.frame(mtx)
y <- seu_test[[class_column]]
y <- ifelse(y == positive_class, positive_class, "ZZZ")
y_test <- factor(y)

data <- new("DataClass", x_train=x_train, y_train=y_train, x_test=x_test, y_test=y_test)
save(data, file = snakemake@output[["data_split"]])

mtx <- seu_train@assays$RNA@data[input_features, ]
mtx <- t(mtx)
x_train_full <- as.data.frame(mtx)
y <- seu_train[[class_column]]
y <- ifelse(y == positive_class, positive_class, "ZZZ")
y_train_full <- factor(y)

### Perform feature selection (scGeneRanger)
scGR_results <- sc_RF_RFE(x=data@x_train, y=data@y_train, folds = fold_indices)
if (scGR_results$optsize > 500) {
  size <- getMinimalSubsetWithinTolerance(scGR_results, tolerance = 0.01)
  selected_features <- getAverageFeatureImps(scGR_results$variables, size)
} else if (scGR_results$optsize > 150) {
  size <- getMinimalSubsetWithinTolerance(scGR_results, tolerance = 0.005)
  selected_features <- getAverageFeatureImps(scGR_results$variables, size)
} else {
  selected_features <- getAverageFeatureImps(scGR_results$variables, scGR_results$optsize)
}

### Save backwards feature elimination cv performance plot
subset_perf_plot <- plot(scGR_results, type = c("g", "o"))
ggsave(file.path(snakemake@output[["cv_results"]]), plot = as.ggplot(subset_perf_plot), width = 6, height = 4, dpi = 300)

### Train final model on full TrainV1V2 with final feature selection
set.seed(42)
xy_train <- cbind(x_train_full, y = y_train_full)
balanced_train_data <- downSample(xy_train, y_train_full, yname = "y")
dim(balanced_train_data)
colnames(balanced_train_data) <- gsub("-", ".", colnames(balanced_train_data))
scGR_results$fit <- ranger(y ~ ., data = balanced_train_data[, c(selected_features$var, "y")], importance = "impurity")

### Assess performance on independent test set
colnames(data@x_test) <- gsub("-", ".", colnames(data@x_test))
predictions <- predict(scGR_results$fit, data = data@x_test[, selected_features$var], type = "response")$predictions
cat(paste("Confusion matrix of final RF model using", length(selected_features$var), "on Test set:", "\n"))
cm <- confusionMatrix(predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_RF <- cm$overall[["Accuracy"]]
F1_RF <- cm$byClass[["F1"]]

### Train kNN classifier for comparison across feature selection methods
custom_knn_with_distance_tie_handling <- function(train, test, cl, k = 5) {
  # Pre-allocate a vector to store predictions
  preds <- vector("character", length = nrow(test))
  
  for (i in seq_len(nrow(test))) {
    # Calculate distances from the i-th test point to all training points
    distances <- sqrt(rowSums((train - test[i, ])^2))
    
    # Find the indices of the k nearest neighbors
    k_indices <- order(distances)[1:k]
    k_distances <- distances[k_indices]
    
    # Handle ties by randomly selecting from equidistant neighbors
    if (length(unique(k_distances)) < k) {
      # Get indices of all points at the cutoff distance
      cutoff_distance <- max(k_distances)
      tie_indices <- which(distances <= cutoff_distance)
      
      # Randomly sample k neighbors among those within the cutoff distance
      selected_indices <- sample(tie_indices, k)
      neighbor_classes <- cl[selected_indices]
    } else {
      # No tie-breaking needed, proceed with nearest k neighbors
      neighbor_classes <- cl[k_indices]
    }
    
    # Determine the predicted class by majority vote among k neighbors
    vote_counts <- table(neighbor_classes)
    max_votes <- max(vote_counts)
    tied_classes <- names(vote_counts)[vote_counts == max_votes]
    
    # Randomly break ties if needed
    preds[i] <- if (length(tied_classes) > 1) {
      sample(tied_classes, 1)
    } else {
      tied_classes
    }
  }
  
  # Return predictions as a factor
  return(factor(preds, levels = unique(cl)))
}


# Define the custom model
knn_custom_model <- list(
  type = "Classification",
  library = "class",
  loop = NULL,
  parameters = data.frame(parameter = "k", class = "numeric", label = "Number of Neighbors"),
  
  fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    list(train = x, cl = y, k = param$k)
  },
  
  predict = function(modelFit, newdata, submodels = NULL) {
    custom_knn_with_distance_tie_handling(train = modelFit$train, test = newdata, cl = modelFit$cl, k = modelFit$k)
  },
  
  prob = NULL,
  
  grid = function(x, y, len = NULL, search = "grid") {
    data.frame(k = seq(5, 23, by = 2))  # Grid for tuning k
  }
)

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

unregister_dopar()

set.seed(42)
train_control <- trainControl(method = "cv", number = 10)
#cl <- makePSOCKcluster(10)
#registerDoParallel(cl)
start_time <- Sys.time()
knn_model <- train(y ~ ., data = balanced_train_data[, c(selected_features$var, "y")], 
                   method = knn_custom_model, trControl = train_control,
                   tuneGrid = data.frame(k = c(5, 7, 9, 11, 13, 15, 17, 19, 23)))
end_time <- Sys.time()
#stopCluster(cl)
print("kNN training finished after:")
print(end_time - start_time)
print(knn_model)

knn_predictions <- predict(knn_model, newdata = data@x_test[, selected_features$var])
cat("Confusion matrix of kNN model on Test set:", "\n")
cm <- confusionMatrix(knn_predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_kNN <- cm$overall[["Accuracy"]]
F1_kNN <- cm$byClass[["F1"]]

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

donor_predictions_rf <- c()
prediction_percentages_rf <- c()
donor_predictions_knn <- c()
prediction_percentages_knn <- c()
for (donor in donors) {
  expr <- FetchData(object = seu_test, vars = donor_column)
  seu_donor <- seu_test[, which(x = expr == donor)]
  if (ncol(seu_donor) < 2) {
    donor_df <- donor_df[donor_df$Donor != donor, ]
    next
  }

  genes_with_dot <- setdiff(gsub(".", "-", selected_features$var, fixed = TRUE), rownames(seu_donor@assays$RNA@data))
  mtx_genes <- intersect(gsub(".", "-", selected_features$var, fixed = TRUE), rownames(seu_donor@assays$RNA@data))
  mtx_genes <- c(mtx_genes, gsub("-", ".", genes_with_dot))
  mtx <- seu_donor@assays$RNA@data[mtx_genes, ]
  mtx <- t(mtx)
  x_donor <- as.data.frame(mtx)
  colnames(x_donor) <- gsub("-", ".", colnames(x_donor))
  y_donor <- seu_donor[[class_column]]
  y_donor <- ifelse(y_donor == positive_class, positive_class, "ZZZ")
  y_donor <- factor(y_donor)

  rf_pred <- predict(scGR_results$fit, data = x_donor, type = "response")$predictions
  pos_pred <- sum(rf_pred == positive_class)
  neg_pred <- sum(rf_pred == "ZZZ")
  if (pos_pred > neg_pred) {
    donor_predictions_rf <- append(donor_predictions_rf, 1)
  }
  else {
    donor_predictions_rf <- append(donor_predictions_rf, 0)
  }
  prediction_percentages_rf <- append(prediction_percentages_rf, pos_pred / (pos_pred + neg_pred))

  knn_pred <- predict(knn_model, newdata = x_donor)
  pos_pred <- sum(knn_pred == positive_class)
  neg_pred <- sum(knn_pred == "ZZZ")
  if (pos_pred > neg_pred) {
    donor_predictions_knn <- append(donor_predictions_knn, 1)
  }
  else {
    donor_predictions_knn <- append(donor_predictions_knn, 0)
  }
  prediction_percentages_knn <- append(prediction_percentages_knn, pos_pred / (pos_pred + neg_pred))
}

donor_df <- cbind(donor_df, donor_predictions_rf, prediction_percentages_rf, donor_predictions_knn, prediction_percentages_knn)
write.csv(donor_df, file = file.path(snakemake@output[["patient_level_predictions"]]), row.names = FALSE)

cat("Patient-level predictions: ", "\n")
cm_donor_pred <- confusionMatrix(factor(donor_predictions_rf, levels = c(1, 0)), factor(donor_df$Group.TrueLabel, levels = c(1, 0)), mode = "everything")
cat("Confusion matrix of RF model patient-level predictions on Test set:", "\n")
print(cm_donor_pred)
Accuracy_RF_Donor_Prediction <- cm_donor_pred$overall[["Accuracy"]]
F1_RF_Donor_Prediction <- cm_donor_pred$byClass[["F1"]]

cm_donor_pred <- confusionMatrix(factor(donor_predictions_knn, levels = c(1, 0)), factor(donor_df$Group.TrueLabel, levels = c(1, 0)), mode = "everything")
cat("Confusion matrix of kNN model patient-level predictions on Test set:", "\n")
print(cm_donor_pred)
Accuracy_kNN_Donor_Prediction <- cm_donor_pred$overall[["Accuracy"]]
F1_kNN_Donor_Prediction <- cm_donor_pred$byClass[["F1"]]

### Save results and model
results_df <- as.data.frame(cbind(selected_features$var, selected_features$Importance))
colnames(results_df) <- c("gene", "imp")
results_df$Accuracy_RF <- Accuracy_RF
results_df$F1_RF <- F1_RF
results_df$Accuracy_RF_Donor_Prediction <- Accuracy_RF_Donor_Prediction
results_df$F1_RF_Donor_Prediction <- F1_RF_Donor_Prediction
results_df$Accuracy_kNN <- Accuracy_kNN
results_df$F1_kNN <- F1_kNN
results_df$Accuracy_kNN_Donor_Prediction <- Accuracy_kNN_Donor_Prediction
results_df$F1_kNN_Donor_Prediction <- F1_kNN_Donor_Prediction

write.csv(results_df, file = snakemake@output[["results"]], row.names = FALSE)
save(scGR_results, file = snakemake@output[["model"]])