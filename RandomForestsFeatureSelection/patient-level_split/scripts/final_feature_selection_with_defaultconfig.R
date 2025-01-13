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
  threshold <- rfe_results$results$WeightedF1[rfe_results$results$Variables == rfe_results$optsize]  - tolerance
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

  cl <- makePSOCKcluster(snakemake@threads)
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

meta_data_column <- snakemake@config[["class_of_interest"]]
donor_slot <- snakemake@config[["donor_slot"]]
group2 <- snakemake@config[["label_of_interest"]]

### Load and subset data
load(snakemake@input[[1]])
seu <- readRDS(snakemake@input[[2]])
seu_train <- subset(seu, cells = c(rownames(data@x_train), rownames(data@x_v2)))
DefaultAssay(seu_train) <- "RNA"
seu_train <- NormalizeData(seu_train)

fold_indices <- CreateGroupFoldIndices(seu_train, sample_col = donor_slot, class_of_interest = meta_data_column, 
                                       label_of_interest = group2, k = 5)
train_rownames <- unique(unlist(lapply(fold_indices, function(x) c(x$train, x$valid))))
seu_subset <- subset(seu_train, cells = train_rownames)
seu_subset
print(sum(seu_subset[[meta_data_column]][, 1] == group2))
print(sum(seu_subset[[meta_data_column]][, 1] != group2))

### Define input feature set
seu_subset <- FindVariableFeatures(seu_subset, selection.method = "vst", nfeatures = 2000)
var_features_top2000 <- VariableFeatures(seu_subset)

### Get expression values and labels
mtx <- seu_subset@assays$RNA@data[var_features_top2000, ]
mtx <- t(mtx)
data@x_train <- as.data.frame(mtx)

y <- seu_subset$Diagnosis
y <- ifelse(y == snakemake@config[["label_of_interest"]], snakemake@config[["label_of_interest"]], "ZZZ")
data@y_train <- factor(y)

gc()

### Feature selection
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

# Test1 ()
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
save(fold_indices, file = snakemake@output[[5]])
