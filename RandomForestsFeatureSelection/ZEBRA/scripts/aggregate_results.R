suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))

set.seed(123)

Acc <- numeric(0)
F1 <- numeric(0)
Acc_RF <- numeric(0)
F1_RF <- numeric(0)

Acc_donor_pred <- numeric(0)
F1_donor_pred <- numeric(0)
Method <- character(0)
Celltype <- character(0)

for (method in c("scGR_defaultconfig", "DeepLift", "FeatureAblation")) {
	for (i in 1:length(snakemake@config[["celltypes"]])) {
		df <- read.csv(snakemake@input[[method]][[i]])
		celltype <- str_extract(snakemake@input[[method]][[i]], "(?<=results/)[^/]+")
		if (method %in% c("DeepLift", "FeatureAblation")) {
			# top 30 features
			Acc <- c(Acc, as.numeric(df$Accuracy_kNN_top30[1]))
			F1 <- c(F1, as.numeric(df$F1_kNN_top30[1]))
			Acc_RF <- c(Acc_RF, as.numeric(df$Accuracy_RF_top30[1]))
			F1_RF <- c(F1_RF, as.numeric(df$F1_RF_top30[1]))
			Method <- c(Method, paste(method, "(top 30)"))
			Celltype <- c(Celltype, celltype)
			# feature set size from scGR with defaultconfig
			Acc <- c(Acc, as.numeric(df$Accuracy_kNN_scGR_defaultconfig[1]))
			F1 <- c(F1, as.numeric(df$F1_kNN_scGR_defaultconfig[1]))
			Acc_RF <- c(Acc_RF, as.numeric(df$Accuracy_RF_scGR_defaultconfig[1]))
			F1_RF <- c(F1_RF, as.numeric(df$F1_RF_scGR_defaultconfig[1]))
			Method <- c(Method, paste(method, "(set size from scGR with default config)"))
			Celltype <- c(Celltype, celltype)
		} else {
			Acc <- c(Acc, as.numeric(df$Accuracy_kNN[1]))
			F1 <- c(F1, as.numeric(df$F1_kNN[1]))
			Acc_RF <- c(Acc_RF, as.numeric(df$Accuracy_RF[1]))
			F1_RF <- c(F1_RF, as.numeric(df$F1_RF[1]))
			Method <- c(Method, method)
			Celltype <- c(Celltype, celltype)
		}
	}
}

results_df <- data.frame(Celltype = Celltype, Method = Method, Acc = Acc, F1 = F1, 
						 Acc_RF = Acc_RF, F1_RF = F1_RF)

results_df <- results_df %>%
  mutate(MethodGroup = case_when(
    grepl("scGR_", Method) ~ "scGeneRanger",
    grepl("DeepLift", Method) ~ "DeepLift",
    grepl("FeatureAblation", Method) ~ "FeatureAblation",
    Method == "DEGs" ~ "DEGs",
    Method == "ElasticNet" ~ "ElasticNet"
  ))

results_df <- results_df %>%
	mutate(MethodLabel = case_when(
		Method == "scGR_defaultconfig" ~ "scGeneRanger (default)",
		Method == "scGR_bestconfig" ~ "scGeneRanger (best)",
		Method == "DEGs" ~ "DEGs",
		Method == "ElasticNet" ~ "ElasticNet",
		Method == "DeepLift (top 30)" ~ "DeepLift (top 30)",
		Method == "DeepLift (set size from scGR with default config)" ~ "DeepLift (*)",
		Method == "DeepLift (set size from scGR with best config)" ~ "DeepLift (**)",
		Method == "FeatureAblation (top 30)" ~ "FeatureAblation (top 30)",
		Method == "FeatureAblation (set size from scGR with default config)" ~ "FeatureAblation (*)",
		Method == "FeatureAblation (set size from scGR with best config)" ~ "FeatureAblation (**)"
	))

write.csv(results_df, file = snakemake@output[[1]])