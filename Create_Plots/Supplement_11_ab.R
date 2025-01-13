library(ggplot2)

source("src/helper.R")
meta_data <- read.csv("data/metadata.csv")
meta_data <- meta_data[meta_data$Visit == 1, ]

colnames(meta_data)

CSF <- c("Quanterix4plex.CSF___NFL","Quanterix4plex.CSF___UCHL1","Quanterix4plex.CSF___Tau",
         "Quanterix4plex.CSF___GFAP","Quanterix3plex.CSF___Ab40","Quanterix3plex.CSF___Ab42" ,"Quanterix3plex.CSF___Tau", 
         "Quanterix.CSF___pTau181")
col_data <- NA
for (meta_col in CSF){
  meta_tmp <- meta_data[, colnames(meta_data) %in% c("SCMD", "Diagnosis", meta_col)]
  colnames(meta_tmp) <- c("SCMD", "Diagnosis", "meta_data")
  meta_tmp$meta <- meta_col
  if (is.na(col_data)[1]) {col_data <- meta_tmp} else {col_data <- rbind(col_data, meta_tmp)}
}
col_data <- col_data[!is.na(col_data$meta_data),]

map <- c("MCI", "PD", "PD-MCI", "AD", "HC")
names(map) <- c("Mild Cognitive Impairment","Parkinson's Disease only","Parkinson's Disease with MCI","Alzheimer's disease" , "Healthy Control")
col_data$Disease_short <- map[col_data$Diagnosis]
col_data$Disease_short <- factor(col_data$Disease_short, levels = c("AD",  "MCI", "PD","PD-MCI","HC"))

plot <- ggplot(col_data, aes(x = Disease_short, y = meta_data, color = Disease_short, fill = Disease_short)) + 
  geom_boxplot() + 
  ggh4x::facet_grid2(.~str_replace(meta, ".*___", ""), scales = "free", independent = "y")+ 
  scale_color_manual(values = darken(color_v))+ scale_fill_manual(values = color_v)+ xlab("") + ylab("")+ theme_adrc()+ 
  labs(fill = "Diagnosis", color = "Diagnosis")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("figures/ComparisonCSF.svg", plot, width = 18, height = 4,unit = "cm")


unique(corr_data$meta[corr_data$p < 0.00002])

brain <- c("Left.vessel","Left.Inf.Lat.Vent","X3rd.Ventricle","Right.Inf.Lat.Vent","Right.vessel","Right.Pallidum")

col_data <- NA
for (meta_col in brain){
  meta_tmp <- meta_data[, colnames(meta_data) %in% c("SCMD", "Diagnosis", meta_col)]
  colnames(meta_tmp) <- c("SCMD", "Diagnosis", "meta_data")
  meta_tmp$meta <- meta_col
  if (is.na(col_data)[1]) {col_data <- meta_tmp} else {col_data <- rbind(col_data, meta_tmp)}
}
col_data <- col_data[!is.na(col_data$meta_data),]


map <- c("MCI", "PD", "PD-MCI", "AD", "HC")
names(map) <- c("Mild Cognitive Impairment","Parkinson's Disease only","Parkinson's Disease with MCI","Alzheimer's disease" , "Healthy Control")
col_data$Disease_short <- map[col_data$Diagnosis]
col_data$Disease_short <- factor(col_data$Disease_short, levels = c("AD",  "MCI", "PD","PD-MCI","HC"))


plot <- ggplot(col_data, aes(x = Disease_short, y = meta_data, color = Disease_short, fill = Disease_short)) + 
  geom_boxplot() + 
  ggh4x::facet_grid2(.~str_replace_all(str_replace(meta, ".*___", ""), "\\.", " "), scales = "free", independent = "y")+ 
  scale_color_manual(values = darken(color_v))+ scale_fill_manual(values = color_v)+
  labs(fill = "Diagnosis", color = "Diagnosis")+
  xlab("") + ylab("")+ theme_adrc()+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("figures/Comparisonbrain.svg", plot, width = 18, height = 4,unit = "cm")
