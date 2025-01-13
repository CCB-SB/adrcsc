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
set.seed(42)
xy_train <- cbind(data@x_train, y = data@y_train)
colnames(xy_train) <- gsub("-", ".", colnames(xy_train))
balanced_train_data <- downSample(xy_train, data@y_train, yname = "y")
dim(balanced_train_data)
rf_model <- ranger(y ~ ., data = balanced_train_data[, c(DEGs_vector, "y")], importance = "impurity")

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

genes_with_dot <- setdiff(gsub(".", "-", DEGs_vector, fixed = TRUE), rownames(seu_test2@assays$RNA@data))
mtx_genes <- intersect(gsub(".", "-", DEGs_vector, fixed = TRUE), rownames(seu_test2@assays$RNA@data))
mtx_genes <- c(mtx_genes, gsub("-", ".", genes_with_dot))
mtx <- seu_test2@assays$RNA@data[mtx_genes, ]
mtx <- t(mtx)
x_test2 <- as.data.frame(mtx)
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
set.seed(123)
train_control <- trainControl(method = "cv", number = 10)
#cl <- makePSOCKcluster(snakemake@threads)
#registerDoParallel(cl)
start_time <- Sys.time()
knn_model <- train(y ~ ., data = balanced_train_data[, c(DEGs_vector, "y")], 
                   method = "knn", trControl = train_control,
                   tuneLength = 10)
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