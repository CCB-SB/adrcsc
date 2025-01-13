suppressMessages(library(Seurat))
suppressMessages(library(rhdf5))
suppressMessages(library(HDF5Array))
suppressMessages(library(dplyr))
suppressMessages(library(caret))

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile, type = "output")
sink(logfile, type = "message")

`%ni%` <- Negate(`%in%`)

seu <- readRDS(file = snakemake@input[[1]])
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu)

load(snakemake@input[[2]])
seu_train <- subset(seu, cells = c(rownames(data@x_train), rownames(data@x_v2)))

expr_matrix <- GetAssayData(seu_train, slot = "data")
seu_train <- SetIdent(seu_train, value = "Diagnosis")

### Filtering of lowly expressed genes
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

# Subset the Seurat object to keep only the pre-filtered genes
seu_train_filtered <- subset(seu_train, features = names(genes_to_keep[genes_to_keep]))

seu_train_filtered



write_h5_DL <- function(exprs_list, h5file_list) {
  
  if (length(unique(lapply(exprs_list, rownames))) != 1) {
    stop("rownames of exprs_list are not identical.")
  }
  
  for (i in seq_along(exprs_list)) {
    if (file.exists(h5file_list[i])) {
      warning("h5file exists! will rewrite it.")
      system(paste("rm", h5file_list[i]))
    }
    
    h5createFile(h5file_list[i])
    h5createGroup(h5file_list[i], "matrix")
    writeHDF5Array(t((exprs_list[[i]])), h5file_list[i], name = "matrix/data")
    h5write(rownames(exprs_list[[i]]), h5file_list[i], name = "matrix/features")
    h5write(colnames(exprs_list[[i]]), h5file_list[i], name = "matrix/barcodes")
    print(h5ls(h5file_list[i]))
  }
}

write_csv_DL <- function(cellType_list, csv_list) {
  
  for (i in seq_along(cellType_list)) {
    if (file.exists(csv_list[i])) {
      warning("csv_list exists! will rewrite it.")
      system(paste("rm", csv_list[i]))
    }
    
    names(cellType_list[[i]]) <- NULL
    write.csv(cellType_list[[i]], file = csv_list[i])
  }
}

downsample_matrix <- function(data_matrix, labels) {
  # Ensure labels are a factor
  labels <- as.factor(labels)
  
  # Identify the minimum class size for downsampling
  min_class_size <- min(table(labels))
  
  # Initialize a list to store indices for each downsampled class
  sampled_indices <- lapply(levels(labels), function(class) {
    # Get indices of the current class
    class_indices <- which(labels == class)
    # Randomly sample without replacement
    sample(class_indices, min_class_size)
  })
  
  # Combine indices from all classes
  downsampled_indices <- unlist(sampled_indices)
  
  # Subset data and labels based on downsampled indices
  balanced_matrix <- data_matrix[downsampled_indices, , drop = FALSE]
  balanced_labels <- labels[downsampled_indices]
  
  # Return matrix and label vector
  list(balanced_matrix = balanced_matrix, balanced_labels = as.vector(balanced_labels))
}

condition_labels <- seu_train_filtered$Diagnosis
train_matrix <- as.matrix(GetAssayData(seu_train_filtered, slot = "data"))

set.seed(42)
balanced_data <- downsample_matrix(t(train_matrix), condition_labels)
balanced_matrix <- t(balanced_data$balanced_matrix)
balanced_labels <- balanced_data$balanced_labels

print(dim(balanced_matrix))
print(table(balanced_labels))

write_h5_DL(exprs_list = list(rna = balanced_matrix), h5file_list = c(snakemake@output[[1]]))
write_csv_DL(cellType_list = list(rna = balanced_labels), csv_list = c(snakemake@output[[2]]))