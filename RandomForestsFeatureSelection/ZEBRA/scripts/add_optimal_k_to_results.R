library(dplyr)

# Define helper function to parse log files for "k"
extract_optimal_k <- function(filepath, num_methods = 1) {
  # Read the log file
  log_lines <- readLines(filepath)
  
  # Find all lines containing the optimal k information
  k_lines <- grep("The final value used for the model was k =", log_lines, value = TRUE)
  
  # Extract the numeric value of k
  k_values <- as.numeric(sub(".*k = (\\d+).*", "\\1", k_lines))
  
  # Handle cases with multiple methods per log file
  if (length(k_values) != num_methods) {
    stop(paste("Unexpected number of k values in", filepath))
  }
  
  return(k_values)
}

# Filepaths from the Snakemake input
input_files <- list(
  scGR_defaultconfig = snakemake@input[["scGR_defaultconfig"]],
  DeepLift = snakemake@input[["DeepLift"]],
  FeatureAblation = snakemake@input[["FeatureAblation"]]
)

# Read in the existing results table
results_df <- read.csv(snakemake@input[["results_df"]])

# Initialize a list to store optimal k values
optimal_k_values <- list()

# Process each method
for (method in names(input_files)) {
  for (file in input_files[[method]]) {
    
    if (method %in% c("DeepLift", "FeatureAblation")) {
      celltype <- basename(dirname(dirname(file)))
      k_values <- extract_optimal_k(file, num_methods = 2)
      
      # Assign specific method names
      if (method == "DeepLift") {
        optimal_k_values <- append(optimal_k_values, list(
          list(Method = "DeepLift (top 30)", Celltype = celltype, k = k_values[1]),
          list(Method = "DeepLift (set size from scGR with default config)", Celltype = celltype, k = k_values[2])
        ))
      } else if (method == "FeatureAblation") {
        optimal_k_values <- append(optimal_k_values, list(
          list(Method = "FeatureAblation (top 30)", Celltype = celltype, k = k_values[1]),
          list(Method = "FeatureAblation (set size from scGR with default config)", Celltype = celltype, k = k_values[2])
        ))
      }
    } else {
      # General case: 1 method per log file
      celltype <- basename(dirname(file))
      k_value <- extract_optimal_k(file, num_methods = 1)
      optimal_k_values <- append(optimal_k_values, list(list(Method = method, Celltype = celltype, k = k_value[1])))
    }
  }
}

# Combine extracted k values into a data frame
optimal_k_df <- do.call(rbind, lapply(optimal_k_values, as.data.frame))

# Join the k values with the results data frame
results_df <- results_df %>%
  left_join(optimal_k_df, by = c("Method", "Celltype")) %>%
  rename(Optimal_k = k)

# Save the updated results to the output file
write.csv(results_df, snakemake@output[[1]], row.names = FALSE)