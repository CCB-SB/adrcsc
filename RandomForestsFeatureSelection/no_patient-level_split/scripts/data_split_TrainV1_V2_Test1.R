library(dplyr)
library(caret)
library(Seurat)

set.seed(123)

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile, type = "output")
sink(logfile, type = "message")

seu <- readRDS(file = snakemake@input[[1]])
DefaultAssay(seu) <- "RNA"
#seu <- NormalizeData(seu)

### Subset to celltype
message("Complete seurat object: ")
seu

annotation_slot <- snakemake@config[["annotation_slot"]]
expr <- FetchData(object = seu, vars = annotation_slot)
seu <- seu[, which(x = expr == snakemake@config[["celltypes"]][[snakemake@params[["celltype"]]]][["name"]])]

message(snakemake@params[["celltype"]], " subset object: ")
seu
message("# ", snakemake@config[["label_of_interest"]], " cells:")
message(sum(seu[[snakemake@config[["class_of_interest"]]]] == snakemake@config[["label_of_interest"]]))
message("# other cells:")
message(sum(seu[[snakemake@config[["class_of_interest"]]]] != snakemake@config[["label_of_interest"]])) 

ExtractData <- function(seu, response_var) {
  data <- seu@assays$RNA@data
  data <- t(data) #%>% as.data.frame() 
  data <- data[, VariableFeatures(seu)]
  gc()
  y <- seu[[response_var]]
  data <- cbind(data, y)
  data[, response_var] <- factor(data[, response_var])
  print(dim(data))
  data
}

### Downsampling function for one-vs-all classification (or binary classification)
DownsampleDataForOneVsAll <- function(data, class_of_interest, label_of_interest, n) {
  labels <- unique(data[, class_of_interest])
  downsampled_data <- data.frame()
  selection_size <- round(n / (length(labels)-1))
  for (label in labels) {
    if (label == label_of_interest) {
      if (sum(data[, class_of_interest] == label_of_interest) > n) {
        selected_indices <- sample(nrow(data[data[, class_of_interest] == label,]), n)
        downsampled_data <- rbind(downsampled_data, data[data[, class_of_interest] == label,][selected_indices,])
      }
      else {
        downsampled_data <- rbind(downsampled_data, data[data[, class_of_interest] == label,])
      }
    }
    else {
      if (sum(data[, class_of_interest] == label) > selection_size) {
        selected_indices <- sample(nrow(data[data[, class_of_interest] == label,]), selection_size)
        downsampled_data <- rbind(downsampled_data, data[data[, class_of_interest] == label,][selected_indices,])
      }
      else {
        downsampled_data <- rbind(downsampled_data, data[data[, class_of_interest] == label,])
      }
    }
  }
  downsampled_data
}

### Helper function that sets binary labels for one-vs-all classification
SetBinaryLabels <- function(data, class_of_interest, label_of_interest) {
  data$binary_label <- "ZZZ"
  data[data[, class_of_interest] == label_of_interest, "binary_label"] <- label_of_interest
  data$binary_label <- factor(data$binary_label)
  data
}

### Custom data class (old)
DataClassOld <- setClass("DataClassOld", slots = c(x_train = "data.frame", y_train = "factor", x_test = "data.frame", y_test = "factor"))

### Wrapper function that prepares the normalized and scaled data from a seurat object for one-vs-all (or binary) classification, combining downsampling and splitting into training and test set
PrepareData <- function(seu, class_of_interest, label_of_interest, n, p = 0.75) {
  
  full_data <- ExtractData(seu, response_var = class_of_interest)
  full_data <- DownsampleDataForOneVsAll(full_data, class_of_interest, label_of_interest, n)
  full_data <- SetBinaryLabels(full_data, class_of_interest, label_of_interest)
  
  x <- select(full_data, -all_of(class_of_interest), -binary_label)
  y <- full_data[, "binary_label"]
  y <- factor(y)
  print(levels(y))
  
  split_indices <- createDataPartition(y, p = p)
  train_indices <- split_indices$Resample1
  test_indices <- setdiff(seq_len(nrow(x)), train_indices)
  
  x_train <- x[train_indices, ]
  x_test <- x[test_indices, ]
  y_train <- y[train_indices]
  y_test <- y[test_indices]
  
  data <- new("DataClassOld", x_train = x_train, y_train = y_train, x_test = x_test, y_test = y_test)
  data
}

### Custom data class (new)
DataClass <- setClass("DataClass", slots = c(x_train = "data.frame", y_train = "factor", x_v2 = "data.frame", y_v2 = "factor", x_test = "data.frame", y_test = "factor"))

data <- PrepareData(seu, class_of_interest = snakemake@config[["class_of_interest"]], label_of_interest = snakemake@config[["label_of_interest"]], 
					n = min(sum(seu[[snakemake@config[["class_of_interest"]]]] == snakemake@config[["label_of_interest"]]), sum(seu[[snakemake@config[["class_of_interest"]]]] != snakemake@config[["label_of_interest"]])))

print("Dimensions of TrainV1 and V2 combined: ")
print(dim(data@x_train))
print("# of AD cells: ")
print(sum(data@y_train == "Alzheimer's disease"))
print("# of CT cells: ")
print(sum(data@y_train == "ZZZ"))

print("Dimensions of Test1: ")
print(dim(data@x_test))
print("# of AD cells: ")
print(sum(data@y_test == "Alzheimer's disease"))
print("# of CT cells: ")
print(sum(data@y_test == "ZZZ"))

split_indices <- createDataPartition(data@y_train, p = 0.666666666666)
trainV1_indices <- split_indices$Resample1
V2_indices <- setdiff(seq_len(nrow(data@x_train)), trainV1_indices)

data <- new("DataClass", x_train = data@x_train[trainV1_indices, ], y_train = data@y_train[trainV1_indices],
			x_v2 = data@x_train[V2_indices, ], y_v2 = data@y_train[V2_indices], x_test = data@x_test, y_test = data@y_test)

print("New TrainV1 / V2 split:")

print("Dimensions of TrainV1: ")
print(dim(data@x_train))
print("# of AD cells: ")
print(sum(data@y_train == "Alzheimer's disease"))
print("# of CT cells: ")
print(sum(data@y_train == "ZZZ"))

print("Dimensions of V2: ")
print(dim(data@x_v2))
print("# of AD cells: ")
print(sum(data@y_v2 == "Alzheimer's disease"))
print("# of CT cells: ")
print(sum(data@y_v2 == "ZZZ"))

save(data, file = snakemake@output[[1]])