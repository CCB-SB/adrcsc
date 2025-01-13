suppressMessages(library(dplyr))
suppressMessages(library(stringr))

set.seed(123)

Acc_test1 <- numeric(0)
F1_test1 <- numeric(0)
Acc_test2 <- numeric(0)
F1_test2 <- numeric(0)
Acc_RF_test1 <- numeric(0)
F1_RF_test1 <- numeric(0)
Acc_RF_test2 <- numeric(0)
F1_RF_test2 <- numeric(0)
Acc_test2_donor_pred <- numeric(0)
F1_test2_donor_pred <- numeric(0)
Method <- character(0)
Celltype <- character(0)

for (method in c("scGR_defaultconfig", "scGR_bestconfig", "DEGs", "ElasticNet", "DeepLift", "FeatureAblation")) {
	for (i in 1:length(snakemake@config[["celltypes"]])) {
		df <- read.csv(snakemake@input[[method]][[i]])
		celltype <- str_extract(snakemake@input[[method]][[i]], "(?<=results/)[^/]+")
		if (method %in% c("DeepLift", "FeatureAblation")) {
			# top 30 features
			Acc_test1 <- c(Acc_test1, as.numeric(df$Accuracy_kNN_Test1_top30[1]))
			F1_test1 <- c(F1_test1, as.numeric(df$F1_kNN_Test1_top30[1]))
			Acc_test2 <- c(Acc_test2, as.numeric(df$Accuracy_kNN_Test2_top30[1]))
			F1_test2 <- c(F1_test2, as.numeric(df$F1_kNN_Test2_top30[1]))
			Acc_RF_test1 <- c(Acc_RF_test1, as.numeric(df$Accuracy_RF_Test1_top30[1]))
			F1_RF_test1 <- c(F1_RF_test1, as.numeric(df$F1_RF_Test1_top30[1]))
			Acc_RF_test2 <- c(Acc_RF_test2, as.numeric(df$Accuracy_RF_Test2_top30[1]))
			F1_RF_test2 <- c(F1_RF_test2, as.numeric(df$F1_RF_Test2_top30[1]))
			Method <- c(Method, paste(method, "(top 30)"))
			Celltype <- c(Celltype, celltype)
			# feature set size from scGR with defaultconfig
			Acc_test1 <- c(Acc_test1, as.numeric(df$Accuracy_kNN_Test1_scGR_defaultconfig[1]))
			F1_test1 <- c(F1_test1, as.numeric(df$F1_kNN_Test1_scGR_defaultconfig[1]))
			Acc_test2 <- c(Acc_test2, as.numeric(df$Accuracy_kNN_Test2_scGR_defaultconfig[1]))
			F1_test2 <- c(F1_test2, as.numeric(df$F1_kNN_Test2_scGR_defaultconfig[1]))
			Acc_RF_test1 <- c(Acc_RF_test1, as.numeric(df$Accuracy_RF_Test1_scGR_defaultconfig[1]))
			F1_RF_test1 <- c(F1_RF_test1, as.numeric(df$F1_RF_Test1_scGR_defaultconfig[1]))
			Acc_RF_test2 <- c(Acc_RF_test2, as.numeric(df$Accuracy_RF_Test2_scGR_defaultconfig[1]))
			F1_RF_test2 <- c(F1_RF_test2, as.numeric(df$F1_RF_Test2_scGR_defaultconfig[1]))
			Method <- c(Method, paste(method, "(set size from scGR with default config)"))
			Celltype <- c(Celltype, celltype)
			# feature set size from scGR with bestconfig
			Acc_test1 <- c(Acc_test1, as.numeric(df$Accuracy_kNN_Test1_scGR_bestconfig[1]))
			F1_test1 <- c(F1_test1, as.numeric(df$F1_kNN_Test1_scGR_bestconfig[1]))
			Acc_test2 <- c(Acc_test2, as.numeric(df$Accuracy_kNN_Test2_scGR_bestconfig[1]))
			F1_test2 <- c(F1_test2, as.numeric(df$F1_kNN_Test2_scGR_bestconfig[1]))
			Acc_RF_test1 <- c(Acc_RF_test1, as.numeric(df$Accuracy_RF_Test1_scGR_bestconfig[1]))
			F1_RF_test1 <- c(F1_RF_test1, as.numeric(df$F1_RF_Test1_scGR_bestconfig[1]))
			Acc_RF_test2 <- c(Acc_RF_test2, as.numeric(df$Accuracy_RF_Test2_scGR_bestconfig[1]))
			F1_RF_test2 <- c(F1_RF_test2, as.numeric(df$F1_RF_Test2_scGR_bestconfig[1]))
			Method <- c(Method, paste(method, "(set size from scGR with best config)"))
			Celltype <- c(Celltype, celltype)
		} else {
			Acc_test1 <- c(Acc_test1, as.numeric(df$Accuracy_kNN_Test1[1]))
			F1_test1 <- c(F1_test1, as.numeric(df$F1_kNN_Test1[1]))
			Acc_test2 <- c(Acc_test2, as.numeric(df$Accuracy_kNN_Test2[1]))
			F1_test2 <- c(F1_test2, as.numeric(df$F1_kNN_Test2[1]))
			Acc_RF_test1 <- c(Acc_RF_test1, as.numeric(df$Accuracy_RF_Test1[1]))
			F1_RF_test1 <- c(F1_RF_test1, as.numeric(df$F1_RF_Test1[1]))
			Acc_RF_test2 <- c(Acc_RF_test2, as.numeric(df$Accuracy_RF_Test2[1]))
			F1_RF_test2 <- c(F1_RF_test2, as.numeric(df$F1_RF_Test2[1]))
			Method <- c(Method, method)
			Celltype <- c(Celltype, celltype)
		}
	}
}

print(Celltype)
print(Method)
print(Acc_test1)
print(F1_test1)
print(Acc_test2)
print(F1_test2)

results_df <- data.frame(Celltype = Celltype, Method = Method, Acc_test1 = Acc_test1, F1_test1 = F1_test1,
						 Acc_test2 = Acc_test2, F1_test2 = F1_test2, Acc_RF_test1 = Acc_RF_test1, F1_RF_test1 = F1_RF_test1, Acc_RF_test2 = Acc_RF_test2, F1_RF_test2 = F1_RF_test2)

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