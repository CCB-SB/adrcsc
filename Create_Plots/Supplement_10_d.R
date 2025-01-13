library(Seurat)
library(data.table)
library(ggpubr)

set.seed(42)

pbmc = readRDS("/local/s8frgran/ADRC/CompleteObjectAnnotated_onlyFirstVisit.rds")
# pbmc2 <- readRDS("/local/s8frgran/ADRC/CompleteObjectAnnotated_secondVisit.rds")
# pbmc3 <- readRDS("/local/s8frgran/ADRC/CompleteObjectAnnotated_thirdVisit.rds")
# 
# pbmc <- merge(pbmc, pbmc2)
# pbmc <- merge(pbmc, pbmc3)
metadata = fread("/local/s8frgran/ADRC/Server/metadata.csv")

cellt_prop <- NA
for (sample in unique(pbmc$PIDN)){
  prop <- as.data.frame((table(pbmc$Manual_Annotation[pbmc$PIDN == sample])))
  prop$Freq <- prop$Freq / sum(prop$Freq)
  data <- data.frame(PIDN = sample)
  for (cellt in unique(pbmc$Manual_Annotation)){
    if(cellt %in% prop$Var1) {data[cellt] <- prop$Freq[prop$Var1 == cellt]} else data[cellt] <- NA
  }
  data$Patient <- pbmc$SCMD[pbmc$PIDN == sample][1]
  data$Visit <- pbmc$Visit[pbmc$PIDN == sample][1]
  if(length(cellt_prop)<2) {cellt_prop <- data} else {cellt_prop <- rbind(cellt_prop, data)}
}

metadata_new <- merge(metadata, cellt_prop, by = c("PIDN"))

write.csv(metadata_new, "SourceData/Supplement_10_d.csv")
#get Celltype-proportions for patients at different Time-Points
df_patients <- data.frame(metadata_new$SCMD, metadata_new$Visit.x)
df_patients <- df_patients[!is.na(metadata_new$`CD14 Monocytes`),]
df_patients <- df_patients[!duplicated(df_patients),]
count <- as.data.frame(table(df_patients$metadata_new.SCMD))
patients_with_moreTP <- as.character(count$Var1[count$Freq>1])

p0 <- metadata_new[metadata_new$SCMD %in% patients_with_moreTP & metadata_new$Visit.x == 1,]
p0 <- p0[!duplicated(p0$SCMD),]
p1 <- metadata_new[metadata_new$SCMD %in% patients_with_moreTP & metadata_new$Visit.x == 2,]
p1 <- p1[!duplicated(p1$SCMD),]

combined_data <- NA
for (celltype in unique(pbmc$Manual_Annotation)){
  df <- (data.frame(x = p0[[celltype]], y = p1[[celltype]], celltype = celltype, type = "timepoints"))
  if (is.na(combined_data)[1]) {combined_data <- df} else {combined_data <- rbind(combined_data, df)}
}


# get Celltype-proportions for different patients from teh same group

patients_with_1TP <- as.character(count$Var1[count$Freq==1])
patients_with_1TP <- unique(metadata_new$SCMD[metadata_new$SCMD %in% patients_with_1TP & metadata_new$Diagnosis == "Healthy Control"])
p0 <- sample(x=patients_with_1TP, size=floor(0.5*length(patients_with_1TP)))
p0 <- metadata_new[metadata_new$SCMD %in% p0,]
p0 <- p0[!duplicated(p0$SCMD),]
p1 <- sample(x= patients_with_1TP[!patients_with_1TP %in% p0], size=floor(0.5*length(patients_with_1TP)))
p1 <- metadata_new[metadata_new$SCMD %in% p1,]
p1 <- p1[!duplicated(p1$SCMD),]


for (celltype in unique(pbmc$Manual_Annotation)){
  df <- (data.frame(x = p0[[celltype]], y = p1[[celltype]], celltype = celltype, type = "same group"))
  if (is.na(combined_data)[1]) {combined_data <- df} else {combined_data <- rbind(combined_data, df)}
}

# get Celltype-proportions for different patients from different groups
patients_with_1TP <- as.character(count$Var1[count$Freq==1])
p0 <- metadata_new[metadata_new$SCMD %in% patients_with_1TP & metadata_new$Diagnosis == "Healthy Control",]
p0 <- p0[!duplicated(p0$SCMD),]
p1 <- metadata_new[metadata_new$SCMD %in% patients_with_1TP & metadata_new$Diagnosis == "Parkinson's Disease only",]
p1 <- p1[!duplicated(p1$SCMD),]

p0 <- p0[p0$SCMD %in% sample(x = p0$SCMD, size = length(p1$PIDN)),]

for (celltype in unique(pbmc$Manual_Annotation)){
  df <- (data.frame(x = p0[[celltype]], y = p1[[celltype]], celltype = celltype, type = "different group"))
  if (is.na(combined_data)[1]) {combined_data <- df} else {combined_data <- rbind(combined_data, df)}
}


patients_with_1TP <- as.character(count$Var1[count$Freq==1])
p0 <- metadata_new[metadata_new$SCMD %in% patients_with_1TP & metadata_new$Diagnosis == "Healthy Control",]
p0 <- p0[!duplicated(p0$SCMD),]
p1 <- metadata_new[metadata_new$SCMD %in% patients_with_1TP & metadata_new$Diagnosis == "Alzheimer's disease",]
p1 <- p1[!duplicated(p1$SCMD),]

p0 <- p0[p0$SCMD %in% sample(x = p0$SCMD, size = length(p1$PIDN)),]

for (celltype in unique(pbmc$Manual_Annotation)){
  df <- (data.frame(x = p0[[celltype]], y = p1[[celltype]], celltype = celltype, type = "different group"))
  if (is.na(combined_data)[1]) {combined_data <- df} else {combined_data <- rbind(combined_data, df)}
}

plot <- ggscatter(combined_data, x = "x", y = "y", color = "type", add = "reg.line", conf.int = TRUE, palette = c("black", "darkgray", "red")) 
ggsave("figures/supplement_10.pdf", plot, width = 5, height = 5)

