suppressMessages(library(Seurat))

set.seed(123)

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile, type = "output")
sink(logfile, type = "message")

### Custom data class
DataClass <- setClass("DataClass", slots = c(x_train = "data.frame", y_train = "factor", x_v2 = "data.frame", y_v2 = "factor", x_test = "data.frame", y_test = "factor"))

seu_trainV1 <- readRDS(file = snakemake@input[[1]])
seu_validV2 <- readRDS(file = snakemake@input[[2]])
seu_test1 <- readRDS(file = snakemake@input[[3]])

mtx <- GetAssayData(object = seu_trainV1, assay = "RNA", slot = "data")
mtx <- as.matrix(mtx[1:5, ])
mtx <- t(mtx)
x_train <- as.data.frame(mtx)
y <- seu_trainV1$Diagnosis
y <- ifelse(y == snakemake@config[["label_of_interest"]], snakemake@config[["label_of_interest"]], "ZZZ")
y <- factor(y)
y_train <- y

mtx <- GetAssayData(object = seu_validV2, assay = "RNA", slot = "data")
mtx <- as.matrix(mtx[1:5, ])
mtx <- t(mtx)
x_v2 <- as.data.frame(mtx)
y <- seu_validV2$Diagnosis
y <- ifelse(y == snakemake@config[["label_of_interest"]], snakemake@config[["label_of_interest"]], "ZZZ")
y <- factor(y)
y_v2 <- y

mtx <- GetAssayData(object = seu_test1, assay = "RNA", slot = "data")
mtx <- as.matrix(mtx[1:5, ])
mtx <- t(mtx)
x_test <- as.data.frame(mtx)
y <- seu_test1$Diagnosis
y <- ifelse(y == snakemake@config[["label_of_interest"]], snakemake@config[["label_of_interest"]], "ZZZ")
y <- factor(y)
y_test <- y

data <- new("DataClass", x_train = x_train, y_train = y_train,
			x_v2 = x_v2, y_v2 = y_v2, x_test = x_test, y_test = y_test)

save(data, file = snakemake@output[[1]])	