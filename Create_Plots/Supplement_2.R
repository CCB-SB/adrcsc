library(Seurat)
library(viridisLite)
library(cowplot)
library(ggpubr)
library(data.table)
library(pbapply)
library(gghalves)
source("scripts/helper.R")

pbmc = readRDS("data/CompleteObjectAnnotated_final.rds")
pbmc$biogroup_short = factor(name2short[as.character(pbmc$Diagnosis)], levels=diagnosis_order)
  

fig_mito_perc_per_diag = ggplot(pbmc@meta.data, aes(x=biogroup_short, y = percent.mt, fill=biogroup_short)) +
  geom_half_violin(data=pbmc@meta.data[pbmc$Sex == "female",], side="l", color=darken(color_v["Female"]), fill=color_v["Female"], alpha = 0.8) + 
  geom_half_boxplot(data=pbmc@meta.data[pbmc$Sex == "female",], side="l", fill=color_v["Female"], width=0.2, alpha = 0.8, outlier.shape = NA) +
  geom_half_violin(data=pbmc@meta.data[pbmc$Sex == "male",], side="r", color=darken(color_v["Male"]), fill=color_v["Male"], alpha = 0.8) + 
  geom_half_boxplot(data=pbmc@meta.data[pbmc$Sex == "male",], side="r", fill=color_v["Male"], width=0.2, alpha = 0.8, outlier.shape = NA) +
  theme_adrc() + ylab("Mitochondrial genes") + xlab("") +
  scale_y_continuous(labels=scales::percent_format(scale = 1), breaks = scales::pretty_breaks(5)) +
  scale_fill_manual(values = color_v) +
  scale_color_manual(values = darken(color_v)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

fig_genes_per_diag = ggplot(pbmc@meta.data, aes(x=biogroup_short, y = nFeature_RNA, fill=biogroup_short)) +
  geom_half_violin(data=pbmc@meta.data[pbmc$Sex == "female",], side="l", color=darken(color_v["Female"]), fill=color_v["Female"], alpha = 0.8) + 
  geom_half_boxplot(data=pbmc@meta.data[pbmc$Sex == "female",], side="l", fill=color_v["Female"], width=0.2, alpha = 0.8, outlier.shape = NA) +
  geom_half_violin(data=pbmc@meta.data[pbmc$Sex == "male",], side="r", color=darken(color_v["Male"]), fill=color_v["Male"], alpha = 0.8) + 
  geom_half_boxplot(data=pbmc@meta.data[pbmc$Sex == "male",], side="r", fill=color_v["Male"], width=0.2, alpha = 0.8, outlier.shape = NA) +
  theme_cowplot(10) + ylab("Genes per cell") + xlab("") +
  scale_y_continuous(labels=c("1k", "2k", "3k", "4k"), breaks = c(1000, 2000, 3000, 4000)) +
  scale_fill_manual(values = color_v) +
  scale_color_manual(values = darken(color_v)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

fig_umis_per_diag = ggplot(pbmc@meta.data, aes(x=biogroup_short, y = nCount_RNA, fill=biogroup_short)) +
  geom_half_violin(data=pbmc@meta.data[pbmc$Sex == "female",], side="l", color=darken(color_v["Female"]), fill=color_v["Female"], alpha = 0.8) + 
  geom_half_boxplot(data=pbmc@meta.data[pbmc$Sex == "female",], side="l", fill=color_v["Female"], width=0.2, alpha = 0.8, outlier.shape = NA) +
  geom_half_violin(data=pbmc@meta.data[pbmc$Sex == "male",], side="r", color=darken(color_v["Male"]), fill=color_v["Male"], alpha = 0.8) + 
  geom_half_boxplot(data=pbmc@meta.data[pbmc$Sex == "male",], side="r", fill=color_v["Male"], width=0.2, alpha = 0.8, outlier.shape = NA) +
  theme_adrc() + ylab("UMIs per cell") + xlab("") +
  scale_y_log10(labels=c("500","1K","2K","5K","10K","25K","50K"), breaks = c(500,1000,2000,5000,10000,25000,50000)) +
  scale_fill_manual(values = color_v) +
  scale_color_manual(values = darken(color_v)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

v1 <- pbmc$Sample
v2 <- pbmc$Sex
cells_per_sample = as.data.table(pbmc@meta.data %>% count(Sample, Sex)) #data.table(Cells=table(pbmc$Sample), Sex = table(pbmc$Sex))
cells_per_sample[, biogroup:=pbmc$biogroup_short[match(Sample, pbmc$Sample)]]


print("Cells per sample")
print(mean(cells_per_sample$Cells.N))
print(sd(cells_per_sample$Cells.N))

fig_cells_per_sample_per_diag = ggplot(cells_per_sample, aes(x=biogroup, y = n, fill=biogroup)) +
  geom_half_violin(data=cells_per_sample[cells_per_sample$Sex == "female",], side="l", color=darken(color_v["Female"]), fill=color_v["Female"], alpha = 0.8) + 
  geom_half_boxplot(data=cells_per_sample[cells_per_sample$Sex == "female",], side="l", fill=color_v["Female"], width=0.2, alpha = 0.8, outlier.shape = NA) +
  geom_half_violin(data=cells_per_sample[cells_per_sample$Sex == "male",], side="r", color=darken(color_v["Male"]), fill=color_v["Male"], alpha = 0.8) + 
  geom_half_boxplot(data=cells_per_sample[cells_per_sample$Sex == "male",], side="r", fill=color_v["Male"], width=0.2, alpha = 0.8, outlier.shape = NA) +
  theme_adrc() + ylab("Cells per sample") + xlab("") +
  scale_y_continuous(labels=c("0", "2K", "4K", "6K", "8K", "10K"), breaks = c(0, 2000, 4000, 6000, 8000, 10000)) +
  scale_fill_manual(values = color_v) +
  scale_color_manual(values = darken(color_v)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

umis_per_sample = as.data.table(pbmc@meta.data[, c("nCount_RNA", "Sample")])[, list(umis=sum(nCount_RNA)), by="Sample"]
umis_per_sample[, biogroup:=pbmc$biogroup_short[match(Sample, pbmc$Sample)]]
umis_per_sample[, Sex:=pbmc$Sex[match(Sample, pbmc$Sample)]]


fig_umis_per_sample_per_diag = ggplot(umis_per_sample, aes(x=biogroup, y = umis, fill=biogroup)) +
  geom_half_violin(data=umis_per_sample[umis_per_sample$Sex == "female",], side="l", color=darken(color_v["Female"]), fill=color_v["Female"], alpha = 0.8) + 
  geom_half_boxplot(data=umis_per_sample[umis_per_sample$Sex == "female",], side="l", fill=color_v["Female"], width=0.2, alpha = 0.8, outlier.shape = NA) +
  geom_half_violin(data=umis_per_sample[umis_per_sample$Sex == "male",], side="r", color=darken(color_v["Male"]), fill=color_v["Male"], alpha = 0.8) + 
  geom_half_boxplot(data=umis_per_sample[umis_per_sample$Sex == "male",], side="r", fill=color_v["Male"], width=0.2, alpha = 0.8, outlier.shape = NA) +
  theme_adrc() + ylab("UMIs per sample") + xlab("") +
  scale_y_continuous(labels=c("0", "5M", "10M", "15M", "20M"), breaks = c(0, 5000000, 10000000, 15000000, 20000000)) +
  scale_fill_manual(values = color_v) +
  scale_color_manual(values = darken(color_v)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

genes_per_sample = pbsapply(unique(pbmc$Sample), function(s){ sum(Matrix::rowSums(pbmc@assays$RNA@counts[,pbmc$Sample == s]) > 0) }, cl=64)
genes_per_sample_df = data.table(Genes=genes_per_sample, Sample=unique(pbmc$Sample))
genes_per_sample_df[, biogroup:=pbmc$biogroup_short[match(Sample, pbmc$Sample)]]
genes_per_sample_df[, Sex:=pbmc$Sex[match(Sample, pbmc$Sample)]]


print("Genes per sample")
print(mean(genes_per_sample_df$Genes))
print(sd(genes_per_sample_df$Genes))



fig_genes_per_sample_per_diag = ggplot(genes_per_sample_df, aes(x=biogroup, y = Genes, fill=biogroup)) +
  geom_half_violin(data=genes_per_sample_df[genes_per_sample_df$Sex == "female",], side="l", color=darken(color_v["Female"]), fill=color_v["Female"], alpha = 0.8) + 
  geom_half_boxplot(data=genes_per_sample_df[genes_per_sample_df$Sex == "female",], side="l", fill=color_v["Female"], width=0.2, alpha = 0.8, outlier.shape = NA) +
  geom_half_violin(data=genes_per_sample_df[genes_per_sample_df$Sex == "male",], side="r", color=darken(color_v["Male"]), fill=color_v["Male"], alpha = 0.8) + 
  geom_half_boxplot(data=genes_per_sample_df[genes_per_sample_df$Sex == "male",], side="r", fill=color_v["Male"], width=0.2, alpha = 0.8, outlier.shape = NA) +
  theme_adrc() + ylab("Genes per sample") + xlab("") +
  scale_y_continuous(labels=c("10K", "15K", "20K", "25K", "30K", "35K"), breaks = c(10000, 15000, 20000, 25000, 30000, 35000)) +
  scale_fill_manual(values = color_v) +
  scale_color_manual(values = darken(color_v)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())


fig_abcde = plot_grid(fig_umis_per_diag, fig_genes_per_diag, fig_mito_perc_per_diag, fig_cells_per_sample_per_diag, fig_umis_per_sample_per_diag, fig_genes_per_sample_per_diag, nrow = 2)
save_plot("figures/Supplement_2_b.svg", fig_abcde, base_height = 100, base_width = 200, unit="mm")


meta_data <- read.csv("/local/s8frgran/ADRC/Server/metadata.csv")
meta_data <- meta_data[,!colnames(meta_data) %in% c("X", "Date.Processed.at.CG")]
meta_data <- meta_data[!duplicated(meta_data),]

p <- ggboxplot(meta_data[meta_data$Visit == 1,], x = "Sex", y = "Age",
               color = "Sex", palette = color_v,
               add = "jitter", short.panel.labs = T,
               facet.by = c("Diagnosis"))+ theme(legend.position = "none")

ggsave("figures/Supplement_2_a.svg", p, width = 7, height = 5)

# Print adjusted p-values for each diagnosis
for (d in unique(meta_data$Diagnosis)){
  p = t.test(meta_data$Age[meta_data$Visit == 1 & meta_data$Sex == "female" & meta_data$Diagnosis == d], meta_data$Age[meta_data$Visit == 1 & meta_data$Sex == "male"& meta_data$Diagnosis == d])
  print(p.adjust(p$p.value, n= 5, method = "BH"))
  print(d)
}


