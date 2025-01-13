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
suppressMessages(library(class))

set.seed(123)

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile, type = "output")
sink(logfile, type = "message")

### Custom data class
DataClass <- setClass("DataClass", slots = c(x_train = "data.frame", y_train = "factor", x_v2 = "data.frame", y_v2 = "factor", x_test = "data.frame", y_test = "factor"))

### scElasticNetFS
### Function to perform feature selection via elastic net logistic regression
sc_ElasticNet_FS <- function(x, y, alpha = 0.5, foldid, class_weights, nfolds = 5) {
  
  x <- as.matrix(x)
  y <- ifelse(y == "ZZZ", 0, 1) #factor?
  
  cl <- makePSOCKcluster(snakemake@threads)
  registerDoParallel(cl)

  start_time <- Sys.time()
  cvfit <- cv.glmnet(x, y, family = "binomial", alpha = alpha, type.measure = "class", nfolds = nfolds, foldid = foldid, weights = class_weights, standardize = TRUE, parallel=TRUE)
  end_time <- Sys.time()

  stopCluster(cl)
  
  elasticnet_time <- end_time - start_time
  print(elasticnet_time)
  
  cvfit
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
    #n <- min(sum(seu[[class_of_interest]][training_indices, ] == label_of_interest), sum(seu[[class_of_interest]][training_indices, ] != label_of_interest))
    #training_indices <- DownsampleIndicesForOneVsAll(training_indices, seu[[class_of_interest]], label_of_interest, n)

    #n <- min(sum(seu[[class_of_interest]][validation_indices, ] == label_of_interest), sum(seu[[class_of_interest]][validation_indices, ] != label_of_interest))
    #validation_indices <- DownsampleIndicesForOneVsAll(validation_indices, seu[[class_of_interest]], label_of_interest, n)

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

### Load and subset data
load(snakemake@input[[1]])
seu <- readRDS(snakemake@input[[2]])
seu_train <- subset(seu, cells = c(rownames(data@x_train), rownames(data@x_v2)))
DefaultAssay(seu_train) <- "RNA"
seu_train <- NormalizeData(seu_train)

#load(snakemake@input[[4]])
#train_rownames <- unique(unlist(lapply(fold_indices, function(x) c(x$train, x$valid))))
#seu_train <- subset(seu_train, cells = train_rownames)
positive_class <- "Alzheimer's disease"

### Determine balanced group k-fold CV fold indices based on donors
fold_indices <- CreateGroupFoldIndices(seu_train, sample_col = "SCMD", class_of_interest = "Diagnosis", 
                                       label_of_interest = "Alzheimer's disease", k = 5)
train_rownames <- unique(unlist(lapply(fold_indices, function(x) c(x$train, x$valid))))
seu_train <- subset(seu_train, cells = train_rownames)

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

y <- seu_train$Diagnosis
y <- ifelse(y == "Alzheimer's disease", "Alzheimer's disease", "ZZZ")
data@y_train <- factor(y)

gc()

foldid <- rep(0, nrow(data@x_train))
# Loop through each fold in `fold_indices` and assign the fold number
for (fold_num in seq_along(fold_indices)) {
  # Match the validation row names to indices in the data
  validation_indices <- which(rownames(data@x_train) %in% fold_indices[[fold_num]]$valid)
  # Assign the fold number to these matched indices in the foldid vector
  foldid[validation_indices] <- fold_num
}
print(head(foldid))
print(table(foldid))

### Determine class weights to adress class imbalance
minority_class <- ifelse(sum(data@y_train == "Alzheimer's disease") >= sum(data@y_train == "ZZZ"), "ZZZ", "Alzheimer's disease")
majority_class <- ifelse(sum(data@y_train == "Alzheimer's disease") < sum(data@y_train == "ZZZ"), "ZZZ", "Alzheimer's disease")
class_weights <- ifelse(y == minority_class, 1 / table(y)[minority_class], 1 / table(y)[majority_class])

### Run Logistic Regression with ElasticNet Regularization (alpha = 0.5)
xy_train <- cbind(data@x_train, y = data@y_train)
colnames(xy_train) <- gsub("-", ".", colnames(xy_train))
balanced_train_data <- xy_train #downSample(xy_train, data@y_train, yname = "y")
dim(balanced_train_data)
data@x_train <- balanced_train_data %>% select(-y)
data@y_train <- balanced_train_data$y

sc_ElasticNet_FS_results <- sc_ElasticNet_FS(data@x_train, data@y_train, foldid = foldid, class_weights = class_weights)
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

#print(setdiff(gsub(".", "-", selected_features$var, fixed = TRUE), rownames(seu_test2@assays$RNA@data)))
genes_with_dot <- setdiff(gsub(".", "-", input_genes, fixed = TRUE), rownames(seu_test2@assays$RNA@data))
mtx_genes <- intersect(gsub(".", "-", input_genes, fixed = TRUE), rownames(seu_test2@assays$RNA@data))
mtx_genes <- c(mtx_genes, gsub("-", ".", genes_with_dot))
mtx <- seu_test2@assays$RNA@data[mtx_genes, ]
#print(setdiff(gsub(".", "-", selected_features$var, fixed = TRUE), rownames(mtx)))
mtx <- t(mtx)
x_test2 <- as.data.frame(mtx)
#print(setdiff(gsub(".", "-", selected_features$var, fixed = TRUE), colnames(x_test2)))
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
set.seed(42)
balanced_train_data <- downSample(xy_train, data@y_train, yname = "y")

rf_model <- ranger(y ~ ., data = balanced_train_data[, c(names(selected_features_ElasticNet), "y")], importance = "impurity")

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

set.seed(123)
print(length(names(selected_features_ElasticNet)))
train_control <- trainControl(method = "cv", number = 10)
#cl <- makePSOCKcluster(snakemake@threads)
#registerDoParallel(cl)
start_time <- Sys.time()
knn_model <- train(y ~ ., data = balanced_train_data[, c(names(selected_features_ElasticNet), "y")], 
                   method = knn_custom_model, trControl = train_control,
                   tuneGrid = data.frame(k = c(5, 7, 9, 11, 13, 15, 17, 19, 23)))
end_time <- Sys.time()
#stopCluster(cl)
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