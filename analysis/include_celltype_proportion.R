library(Seurat)
library(data.table)

pbmc = readRDS("results/annotated_celltypes/CompleteObjectAnnotated.rds")
metadata = fread("results/metadata_with_new_biomarker_data.csv")

cellt_prop <- NA
for (sample in unique(pbmc$Sample)){
  prop <- as.data.frame((table(pbmc$celltype_cluster[pbmc$Sample == sample])))
  prop$Freq <- prop$Freq / sum(prop$Freq)
  data <- data.frame(Sample = sample)
  for (cellt in unique(pbmc$celltype_cluster)){
    if(cellt %in% prop$Var1) {data[cellt] <- prop$Freq[prop$Var1 == cellt]} else data[cellt] <- NA
  }
  data$Patient <- pbmc$Patient[pbmc$Sample == sample][1]
  data$Visit <- pbmc$Visit[pbmc$Sample == sample][1]
  if(length(cellt_prop)<2) {cellt_prop <- data} else {cellt_prop <- rbind(cellt_prop, data)}
}

metadata_new <- merge(metadata, cellt_prop, by = c("Sample"))
write.csv(metadata_new, "results/metadata_with_celltype_cluster_proportions_and_new_biomarker_data.csv")
