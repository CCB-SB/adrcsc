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

### Load DEGs table and get sig. DEGs for current celltype
DEGs <- read.csv(file = snakemake@input[[1]])
celltype <- snakemake@config[["celltypes"]][[snakemake@params[[1]]]][["name"]]
DEGs <- DEGs %>% filter(cluster_id == celltype, category != "Not deregulated")
DEGs_vector <- gsub("-", ".", DEGs$gene)
cat(celltype, "\n")
cat("Number of significant DEGs: ", length(DEGs_vector), "\n")

### Load and subset data
load(snakemake@input[[2]])
seu <- readRDS(snakemake@input[[3]])
seu_train <- subset(seu, cells = c(rownames(data@x_train), rownames(data@x_v2)))
DefaultAssay(seu_train) <- "RNA"
seu_train <- NormalizeData(seu_train)

### Get expression values and labels
mtx <- seu_train@assays$RNA@data[DEGs$gene, ]
mtx <- t(mtx)
data@x_train <- as.data.frame(mtx)

common_levels <- union(levels(data@y_train), levels(data@y_v2))
y_trainV1V2 <- factor(c(as.character(data@y_train), as.character(data@y_v2)), levels = common_levels)
data@y_train <- y_trainV1V2

gc()

### Train RF model on full TrainV1V2 with DEGs
xy_train <- cbind(data@x_train, y = data@y_train)
colnames(xy_train) <- gsub("-", ".", colnames(xy_train))
rf_model <- ranger(y ~ ., data = xy_train[, c(DEGs_vector, "y")], importance = "impurity")

### Assess performance on test sets Test1 and Test2
seu_test1 <- subset(seu, cells = rownames(data@x_test))
seu_test2 <- readRDS(file = snakemake@input[[4]])
DefaultAssay(seu_test1) <- "RNA"
DefaultAssay(seu_test2) <- "RNA"
seu_test1 <- NormalizeData(seu_test1)
seu_test2 <- NormalizeData(seu_test2)

# Test1 (Test cells from same patients)
genes_with_dot <- setdiff(gsub(".", "-", DEGs_vector, fixed = TRUE), rownames(seu_test1@assays$RNA@data))
mtx_genes <- intersect(gsub(".", "-", DEGs_vector, fixed = TRUE), rownames(seu_test1@assays$RNA@data))
mtx_genes <- c(mtx_genes, gsub("-", ".", genes_with_dot))
mtx <- seu_test1@assays$RNA@data[mtx_genes, ]
mtx <- t(mtx)
data@x_test <- as.data.frame(mtx)
colnames(data@x_test) <- gsub("-", ".", colnames(data@x_test))
gc()

predictions <- predict(rf_model, data = data@x_test, type = "response")$predictions
cat("Confusion matrix of DEGs RF model on Test1:", "\n")
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

#print(setdiff(gsub(".", "-", selected_features$var, fixed = TRUE), rownames(seu_test2@assays$RNA@data)))
genes_with_dot <- setdiff(gsub(".", "-", DEGs_vector, fixed = TRUE), rownames(seu_test2@assays$RNA@data))
mtx_genes <- intersect(gsub(".", "-", DEGs_vector, fixed = TRUE), rownames(seu_test2@assays$RNA@data))
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

predictions <- predict(rf_model, data = x_test2, type = "response")$predictions
cat("Confusion matrix of DEGs RF model on Test2:", "\n")
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
train_control <- trainControl(method = "cv", number = 10)
#cl <- makePSOCKcluster(snakemake@threads)
#registerDoParallel(cl)
start_time <- Sys.time()
knn_model <- train(y ~ ., data = xy_train[, c(DEGs_vector, "y")], 
                   method = knn_custom_model, trControl = train_control,
                   tuneGrid = data.frame(k = c(5, 7, 9, 11, 13, 15, 17, 19, 23)))
end_time <- Sys.time()
#stopCluster(cl)
print("kNN training finished after:")
print(end_time - start_time)
print(knn_model)

knn_predictions <- predict(knn_model, newdata = data@x_test)
cat("Confusion matrix of DEGs kNN model on Test1:", "\n")
cm <- confusionMatrix(knn_predictions, data@y_test, mode = "everything")  
print(cm)
Accuracy_kNN_Test1 <- cm$overall[["Accuracy"]]
F1_kNN_Test1 <- cm$byClass[["F1"]]

knn_predictions <- predict(knn_model, newdata = x_test2)
cat("Confusion matrix of DEGs kNN model on Test2:", "\n")
cm <- confusionMatrix(knn_predictions, y_test2, mode = "everything")  
print(cm)
Accuracy_kNN_Test2 <- cm$overall[["Accuracy"]]
F1_kNN_Test2 <- cm$byClass[["F1"]]

### Save results
DEGs$Accuracy_RF_Test1 <- Accuracy_RF_Test1
DEGs$Accuracy_RF_Test2 <- Accuracy_RF_Test2
DEGs$Accuracy_kNN_Test1 <- Accuracy_kNN_Test1
DEGs$Accuracy_kNN_Test2 <- Accuracy_kNN_Test2
DEGs$F1_RF_Test1 <- F1_RF_Test1
DEGs$F1_RF_Test2 <- F1_RF_Test2
DEGs$F1_kNN_Test1 <- F1_kNN_Test1
DEGs$F1_kNN_Test2 <- F1_kNN_Test2

write.csv(DEGs, file = snakemake@output[[1]], row.names = FALSE)