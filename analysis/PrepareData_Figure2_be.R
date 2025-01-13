library(Seurat)


meta_data <- read.csv("/local/s8frgran/ADRC/Server/metadata.csv")
pbmc <- readRDS("/local/s8frgran/ADRC/CompleteObjectAnnotated_onlyFirstVisit.rds")
#length(unique(pbmc$Sample[pbmc$Sample.ID..From.Stanford. %in% meta_data$Sample.ID..From.Stanford.[meta_data$Visit == 1]]))
pbmc <- subset(pbmc, subset = Sample.ID..From.Stanford. %in% meta_data$Sample.ID..From.Stanford.[meta_data$Visit == 1])

cell_counts_all <- pbmc@meta.data %>% group_by(celltype_cluster, SCMD) %>%  summarise(n = n())

diagn <- c("Mild Cognitive Impairment","Alzheimer's disease","Parkinson's Disease with MCI","Parkinson's Disease only")
names(diagn) <- c("MCI", "AD", "PDMCI", "PD")

meta <- as.data.frame(pbmc@meta.data %>% group_by(SCMD, Sex, Age, Diagnosis, ApoE))
meta <- meta[, colnames(meta) %in% c("SCMD", "Sex", "Age", "Diagnosis", "ApoE")]
meta <- meta[!duplicated(meta),]

meta$Diagnosis <- diagn[meta$Diagnosis]

for (ct in unique(cell_counts_all$celltype_cluster)){
  cell_counts <- cell_counts_all[cell_counts_all$celltype_cluster == ct,]
  cell_counts<- cell_counts[colnames(cell_counts) %in% c("SCMD", "n")]
  print(cell_counts)
  meta <- merge(meta, cell_counts, by = "SCMD", all.x = T)
  colnames(meta)[length(colnames(meta))] <- ct
  meta[[ct]][is.na(meta[[ct]])] <- 0
}

# 
# #meta$Diagnosis[meta$Diagnosis == "Alzheimer's disease"] <- "Z.Disase.AD"
# write.csv(meta, "CellCounts_fine.csv", row.names = F)

meta$ApoE[meta$ApoE == ""] <- "unknown"
write.csv(meta[meta$ApoE != "",], "CellCounts_Apoe_fine.csv", row.names = F)

################################################################################

cell_counts_all <- pbmc@meta.data %>% group_by(celltype, SCMD) %>%  summarise(n = n())

diagn <- c("Mild Cognitive Impairment","Alzheimer's disease","Parkinson's Disease with MCI","Parkinson's Disease only")
names(diagn) <- c("MCI", "AD", "PDMCI", "PD")

meta <- as.data.frame(pbmc@meta.data %>% group_by(SCMD, Sex, Age, Diagnosis, ApoE))
meta <- meta[, colnames(meta) %in% c("SCMD", "Sex", "Age", "Diagnosis", "ApoE")]
meta <- meta[!duplicated(meta),]

meta$Diagnosis <- diagn[meta$Diagnosis]

for (ct in unique(cell_counts_all$celltype)){
  cell_counts <- cell_counts_all[cell_counts_all$celltype == ct,]
  cell_counts<- cell_counts[colnames(cell_counts) %in% c("SCMD", "n")]
  print(cell_counts)
  meta <- merge(meta, cell_counts, by = "SCMD", all.x = T)
  colnames(meta)[length(colnames(meta))] <- ct
  meta[[ct]][is.na(meta[[ct]])] <- 0
}

# 
# #meta$Diagnosis[meta$Diagnosis == "Alzheimer's disease"] <- "Z.Disase.AD"
# write.csv(meta, "CellCounts_fine.csv", row.names = F)

meta$ApoE[meta$ApoE == ""] <- "unknown"
write.csv(meta[meta$ApoE != "",], "CellCounts_Apoe_broad.csv", row.names = F)

