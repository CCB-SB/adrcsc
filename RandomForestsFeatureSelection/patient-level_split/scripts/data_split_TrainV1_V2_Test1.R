suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(caret))

set.seed(123)

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile, type = "output")
sink(logfile, type = "message")

### Load data and fetch config values
seu <- readRDS(file = snakemake@input[[1]])
donor_slot <- snakemake@config[["donor_slot"]]
annotation_slot <- snakemake@config[["annotation_slot"]]
metadata_column <- snakemake@config[["class_of_interest"]]
celltype <- snakemake@config[["celltypes"]][[snakemake@params[[1]]]][["name"]]

### Split donors in TrainV1, V2 and Test1 set
group1_donors <- unique(seu[[donor_slot]][, 1][seu[[metadata_column]][, 1] != snakemake@config[["label_of_interest"]]])
group2_donors <- unique(seu[[donor_slot]][, 1][seu[[metadata_column]][, 1] == snakemake@config[["label_of_interest"]]])

train_donors_group1 <- sample(group1_donors, size = floor(0.5 * length(group1_donors)))
remaining_donors_group1 <- setdiff(group1_donors, train_donors_group1)
validation_donors_group1 <- sample(remaining_donors_group1, size = floor(0.5 * length(remaining_donors_group1)))
test_donors_group1 <- setdiff(remaining_donors_group1, validation_donors_group1)

train_donors_group2 <- sample(group2_donors, size = floor(0.5 * length(group2_donors)))
remaining_donors_group2 <- setdiff(group2_donors, train_donors_group2)
validation_donors_group2 <- sample(remaining_donors_group2, size = floor(0.5 * length(remaining_donors_group2)))
test_donors_group2 <- setdiff(remaining_donors_group2, validation_donors_group2)

train_donors <- c(train_donors_group1, train_donors_group2)
validation_donors <- c(validation_donors_group1, validation_donors_group2)
test_donors <- c(test_donors_group1, test_donors_group2)

print("TrainV1 donors: ")
print(train_donors)
print("V2 donors: ")
print(validation_donors)
print("Test1 donors: ")
print(test_donors)

### Subset to celltype and split data
expr <- FetchData(object = seu, vars = annotation_slot)
seu_celltype <- seu[, which(x = expr == celltype)]
seu_train <- seu_celltype[, seu_celltype[[donor_slot]][, 1] %in% train_donors]
seu_valid <- seu_celltype[, seu_celltype[[donor_slot]][, 1] %in% validation_donors]
seu_test <- seu_celltype[, seu_celltype[[donor_slot]][, 1] %in% test_donors]

print(celltype)
print("Complete Seurat object: ")
seu_celltype
print("Training/V1 set Seurat object: ")
seu_train
print("# of AD cells:")
print(sum(seu_train$Diagnosis == "Alzheimer's disease"))
print("Validation set V2 Seurat object: ")
seu_valid
print("# of AD cells:")
print(sum(seu_valid$Diagnosis == "Alzheimer's disease"))
print("Test1 set Seurat object: ")
seu_test
print("# of AD cells:")
print(sum(seu_test$Diagnosis == "Alzheimer's disease"))

saveRDS(seu_train, file = snakemake@output[[1]])
saveRDS(seu_valid, file = snakemake@output[[2]])
saveRDS(seu_test, file = snakemake@output[[3]])