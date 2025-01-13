colors = fread("colors.csv", strip.white = F)
color_v = colors$Color
names(color_v) = colors$ID



meta_data <- read.csv("/local/s8frgran/ADRC/Server/metadata.csv")
pbmc <- readRDS("/local/s8frgran/ADRC/Upload/Pipeline_All/results_TP1/annotated_celltypes/CompleteObjectAnnotated.rds")
#length(unique(pbmc$Sample[pbmc$Sample.ID..From.Stanford. %in% meta_data$Sample.ID..From.Stanford.[meta_data$Visit == 1]]))
pbmc$Visit <- 1
pbmc$Visit[pbmc$Sample.ID..From.Stanford. %in% meta_data$Sample.ID..From.Stanford.[meta_data$Visit ==2]] <- 2
pbmc$Visit[pbmc$Sample.ID..From.Stanford. %in% meta_data$Sample.ID..From.Stanford.[meta_data$Visit ==3]] <- 3

pbmc <- subset(pbmc, subset = SCMD %in% pbmc$SCMD[pbmc$Visit != 1])

cell_counts_all <- pbmc@meta.data %>% group_by(L2, SCMD, Visit) %>%  summarise(n = n())

meta <- as.data.frame(pbmc@meta.data %>% group_by(SCMD, Sex, Age, Diagnosis, ApoE, Visit))
meta <- meta[, colnames(meta) %in% c("SCMD", "Sex", "Age", "Diagnosis", "ApoE", "Visit")]
meta <- meta[!duplicated(meta),]

for (ct in unique(cell_counts_all$L2)){
  cell_counts <- cell_counts_all[cell_counts_all$L2 == ct,]
  cell_counts<- cell_counts[colnames(cell_counts) %in% c("SCMD", "n", "Visit")]
  print(cell_counts)
  meta <- merge(meta, cell_counts, by = c("SCMD", "Visit"), all.x = T)
  colnames(meta)[length(colnames(meta))] <- ct
  meta[[ct]][is.na(meta[[ct]])] <- 0
}

corr_data <- NA
for (d in unique(meta$Diagnosis)){
  for (ct in unique(cell_counts_all$L2)){
    corr <- cor.test(meta[[ct]][meta$Diagnosis == d], meta$Visit[meta$Diagnosis == d])
    
    data <- data.frame(cor = corr$estimate, p = corr$p.value, Diagnosis = d, celltype = ct )
    if (is.na(corr_data)[1]){corr_data <- data} else {corr_data <- rbind(corr_data, data)}
  }
}

corr_data$celltype <- factor(corr_data$celltype, levels = CellCluster_order)

map <- c("MCI", "PD", "PD-MCI", "AD", "HC")
names(map) <- c("Mild Cognitive Impairment","Parkinson's Disease only","Parkinson's Disease with MCI","Alzheimer's disease" , "Healthy Control")
corr_data$Disease_short <- map[corr_data$Diagnosis]
corr_data$Disease_short <- factor(corr_data$Disease_short, levels = c("AD",  "MCI", "PD","PD-MCI","HC"))

corr_data$p.adj <- p.adjust(corr_data$p)
plot <- ggplot(corr_data, aes(x = celltype, y = Disease_short, fill = cor))+ geom_tile()+ 
  xlab("") + 
  ylab("") + theme_adrc()+
  scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027", midpoint = 0, name = "Pearson's correlation")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


labels= ggplot(corr_data, aes(x=celltype, y=1, fill=celltype)) + geom_tile() +
  scale_fill_manual(values = color_v) +
  theme_void() + guides(fill="none")

plot <- (plot %>% insert_top(labels, height=.25))

ggsave("figures/Supplement_10_e.svg", plot, width = 10, height = 5, unit = "cm")

write.csv(corr_data, "SourceData/Figure_10_e.csv")
