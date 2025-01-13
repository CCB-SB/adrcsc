library(stringr)
library(dplyr)

ntree_values <- c()
maxdepth_values <- c()
feature_set_values <- c()
f1_score_values <- c()
celltype_values <- c()

additional_feature_set_names = c("top6000HVGs", "top8000HVGs", "LowExpressionFilter", "DEGs")
feature_set_names = c("random2000_1", "random2000_2", "random2000_3", "top2000HVGs", "top4000HVGs")

for (i in 1:length(snakemake@config[["celltypes"]])) {
	celltype <- str_extract(snakemake@input[["hyperparams1"]][[i]], "(?<=results/)[^/]+")
	logfile <- readLines(snakemake@input[["hyperparams1"]][[i]])
	for (j in seq_along(logfile)) {
	  # Extract hyperparameter configuration
	  if (str_detect(logfile[j], "Running scGeneRanger")) {
	  	#print(logfile[j])
	    ntree <- as.numeric(str_extract(logfile[j], "(?<=ntree =  )\\d+"))
	    maxdepth <- as.numeric(str_extract(logfile[j], "(?<=maxdepth =  )\\d+"))
	    feature_set <- feature_set_names[as.numeric(str_extract(logfile[j], "(?<=input feature set =  )\\d+"))]
	    
	    # Store the extracted values
	    ntree_values <- c(ntree_values, ntree)
	    maxdepth_values <- c(maxdepth_values, maxdepth)
	    feature_set_values <- c(feature_set_values, feature_set)
	  }
	  
	  # Extract F1 score from confusion matrix
	  if (str_detect(logfile[j], "F1 :")) {
	    f1_score <- as.numeric(str_extract(logfile[j], "(?<=F1 : )\\d+\\.\\d+"))
	    f1_score_values <- c(f1_score_values, f1_score)
	    celltype_values <- c(celltype_values, celltype)
	  }
	}
	logfile <- readLines(snakemake@input[["hyperparams2"]][[i]])
	for (j in seq_along(logfile)) {
	  # Extract hyperparameter configuration
	  if (str_detect(logfile[j], "Running scGeneRanger")) {
	    ntree <- as.numeric(str_extract(logfile[j], "(?<=ntree =  )\\d+"))
	    maxdepth <- as.numeric(str_extract(logfile[j], "(?<=maxdepth =  )\\d+"))
	    feature_set <- additional_feature_set_names[as.numeric(str_extract(logfile[j], "(?<=input feature set =  )\\d+"))]
	    
	    # Store the extracted values
	    ntree_values <- c(ntree_values, ntree)
	    maxdepth_values <- c(maxdepth_values, maxdepth)
	    feature_set_values <- c(feature_set_values, feature_set)
	  }
	  
	  # Extract F1 score from confusion matrix
	  if (str_detect(logfile[j], "F1 :")) {
	    f1_score <- as.numeric(str_extract(logfile[j], "(?<=F1 : )\\d+\\.\\d+"))
	    f1_score_values <- c(f1_score_values, f1_score)
	    celltype_values <- c(celltype_values, celltype)
	  }
	}
} 

# Combine all values into a data frame
results <- data.frame(
  ntree = ntree_values,
  maxdepth = maxdepth_values,
  feature_set = feature_set_values,
  F1_on_V2 = f1_score_values,
  celltype = celltype_values
)

write.csv(results, snakemake@output[[2]])

best_config <- results %>% 
  group_by(ntree, maxdepth, feature_set) %>%
  summarize(mean_F1_on_V2 = mean(F1_on_V2)) %>%
  arrange(desc(mean_F1_on_V2))

write.csv(best_config, snakemake@output[[1]])
