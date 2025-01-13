library(stringr)
library(edgeR)
library(ggplot2)
library(biomaRt)
library(data.table)

md1 <- read.csv("ROSMAP_bulk/ROSMAP_assay_rnaSeq_metadata.csv")
md2 <- read.csv("ROSMAP_bulk/ROSMAP_biospecimen_metadata.csv")
md3 <- read.csv("ROSMAP_bulk/ROSMAP_clinical.csv")

md2 <- md2[md2$specimenID %in% md1$specimenID & md2$assay == "rnaSeq",]
md3 <- md3[md3$individualID %in% md2$individualID,]

meta_data <- merge(md1, md2, by = "specimenID")
meta_data <- merge(meta_data, md3, by = "individualID") 

meta_data <- meta_data[meta_data$organ == "brain",]

length(unique(meta_data$individualID[meta_data$cogdx %in% c(2,3)]))# MCI
length(unique(meta_data$individualID[meta_data$cogdx %in% c(4,5)])) # AD
length(unique(meta_data$individualID[meta_data$cogdx %in% c(1)])) # HC

length((meta_data$individualID[meta_data$cogdx %in% c(2,3)]))# MCI
length((meta_data$individualID[meta_data$cogdx %in% c(4,5)])) # AD
length((meta_data$individualID[meta_data$cogdx %in% c(1)])) # HC

meta_data$individualID[duplicated(meta_data$individualID) & meta_data$cogdx %in% c(2,3)]

length((meta_data$individualID[meta_data$cogdx %in% c(1)& meta_data$tissue == "Head of caudate nucleus"] )) 


gene_expression <- read.csv("ROSMAP_bulk/ROSMAP_RNAseq_FPKM_gene.tsv", sep = "\t")
rownames(gene_expression) <- gene_expression$gene_id
gene_expression <- gene_expression[,-c(1,2)]
sample_names <- colnames(gene_expression) # [-c(1,2)]

sampels <- substr(sample_names, 2, 100)
sampels <- str_replace(sampels, "_[0-9]*$", "")

sample_mapping <- data.frame(sample_names, sampels)

sampels = sampels[sampels %in% meta_data$specimenID]

meta_data = meta_data[meta_data$specimenID %in% sampels,]
sample_mapping = sample_mapping[sample_mapping$sampels %in% meta_data$specimenID,]
sample_mapping = sample_mapping[!duplicated(sample_mapping$sampels),]
gene_expression <- gene_expression[,colnames(gene_expression) %in% sample_mapping$sample_names]
sample_mapping$specimenID <- sample_mapping$sampels

meta_data <- merge(meta_data, sample_mapping, by = "specimenID")
d0 <- DGEList(gene_expression)
# filter ow-expressed genes
cutoff <- 10
drop <- which(apply((gene_expression), 1, max) < cutoff)
d <- d0[-drop,] 


meta_data$Diagnosis <- "Other"
meta_data$Diagnosis[meta_data$cogdx %in% c(4,5)] <- "AD"
meta_data$Diagnosis[meta_data$cogdx %in% c(2,3)] <- "MCI"
meta_data$Diagnosis[meta_data$cogdx %in% c(1)] <- "HC"

diag <- meta_data$Diagnosis[order(match(meta_data$sample_names, colnames(gene_expression)))]
sex <- meta_data$msex[order(match(meta_data$sample_names, colnames(gene_expression)))]
group <- interaction(diag, sex)

mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T, normalize.method = "none")
fit <- lmFit(y, mm)
contr <- makeContrasts(groupAD.0 - groupHC.0, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_AD_fem <- topTable(tmp, sort.by = "P", n = Inf)

plotSA(tmp, main="Final model: Mean-variance trend")

top.table_AD_fem[top.table_AD_fem$adj.P.Val < 0.05,]

top.table_AD_fem$gene_id <- rownames(top.table_AD_fem)
write.csv(top.table_AD_fem, "ROSMAP_bulk/DEGs_female.csv")

contr <- makeContrasts(groupAD.1 - groupHC.1, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_AD_male <- topTable(tmp, sort.by = "P", n = Inf)

plotSA(tmp, main="Final model: Mean-variance trend")

top.table_AD_male[top.table_AD_male$adj.P.Val < 0.05,]
top.table_AD_male$gene_id <- rownames(top.table_AD_male)
write.csv(top.table_AD_male, "ROSMAP_bulk/DEGs_male.csv")

top.table_AD_fem <- read.csv("ROSMAP_bulk/DEGs_female.csv")
top.table_AD_male <- read.csv("ROSMAP_bulk/DEGs_male.csv")

degs_female <- read.csv("Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_L2_with_additional_stats.p_filtered.csv", sep = "\t")
genes_sc <- unique(degs_female$gene_id[degs_female$contrast == "ADvsHC"])
genes_common <- genes_sc[genes_sc %in% rownames(top.table_AD_fem)]

degs_female_subset <- degs_female[degs_female$gene_id %in% genes_common,]
degs_female_subset <- degs_female_subset[,colnames(degs_female_subset) %in% c("gene_id", "cluster_id", "logFC")]
#top.table_AD_fem$gene_id <- rownames(top.table_AD_fem)
top.table_AD_fem <- top.table_AD_fem[top.table_AD_fem$gene_id %in% genes_common,]
top.table_AD_fem$cluster_id <- "bulk Brain"
top.table_AD_fem <- top.table_AD_fem[,colnames(top.table_AD_fem) %in% c("gene_id", "cluster_id", "logFC")]
degs_to_plot <- rbind(degs_female_subset, top.table_AD_fem)
degs_to_plot$logFC[degs_to_plot$cluster_id == "bulk Brain"] <- degs_to_plot$logFC[degs_to_plot$cluster_id == "bulk Brain"] / max(abs(degs_to_plot$logFC[degs_to_plot$cluster_id == "bulk Brain"]))
degs_to_plot$logFC[degs_to_plot$cluster_id != "bulk Brain"] <- degs_to_plot$logFC[degs_to_plot$cluster_id != "bulk Brain"] / max(abs(degs_to_plot$logFC[degs_to_plot$cluster_id != "bulk Brain"]))
ggplot(degs_to_plot, aes(x = gene_id, y = cluster_id, color = logFC))+ geom_point()+ scale_color_gradient2(low="#4575b4", mid="white", high="#d73027")+ theme_classic()



degs_male <- read.csv("Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_L2_with_additional_stats.p_filtered.csv", sep = "\t")
genes_sc <- unique(degs_male$gene_id[degs_male$contrast == "ADvsHC"])
genes_common <- genes_sc[genes_sc %in% (top.table_AD_male$gene_id)]

degs_male_subset <- degs_male[degs_male$gene_id %in% genes_common,]
degs_male_subset <- degs_male_subset[,colnames(degs_male_subset) %in% c("gene_id", "cluster_id", "logFC")]
#top.table_AD_fem$gene_id <- rownames(top.table_AD_fem)
top.table_AD_male <- top.table_AD_male[top.table_AD_male$gene_id %in% genes_common,]
top.table_AD_male$cluster_id <- "bulk Brain"
top.table_AD_male <- top.table_AD_male[,colnames(top.table_AD_male) %in% c("gene_id", "cluster_id", "logFC")]
degs_to_plot <- rbind(degs_male_subset, top.table_AD_male)
ggplot(degs_to_plot, aes(x = gene_id, y = cluster_id, color = logFC))+ geom_point()+ scale_color_gradient2(low="#4575b4", mid="white", high="#d73027", limits = c(-1,1))+ theme_classic()


# ENSG00000218227.3








gene <-  "ENSG00000105643.3" #
gene_expr <- gene_expression[gene,]
gene_expr_meta <- meta_data[order(match(meta_data$sample_names,colnames(gene_expr))),]
tmp <- gene_expr[,order(match(colnames(gene_expr),meta_data$sample_names))]
gene_expr_meta$expr <- unlist(tmp)


ggplot(gene_expr_meta[gene_expr_meta$Diagnosis %in% c("AD", "HC"),], aes(x = Diagnosis, y = expr))+ facet_wrap(~msex)+ geom_boxplot()

httr::set_config(httr::config(ssl_verifypeer=FALSE))
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

ensemblid2gene = fread("Pipeline/data/annotations/geneid2symbol.final.csv")
top.table_AD_fem$gene <- rownames(top.table_AD_fem)
top.table_AD_fem$gene_id = ensemblid2gene$GeneSymbol[match(top.table_AD_fem$gene, ensemblid2gene$Geneid)]

id2biotype = data.table(getBM(c("ensembl_gene_id_version", "transcript_biotype"), filters="ensembl_gene_id_version", values=unique(top.table_AD_fem$gene), mart=ensembl))
id2biotype = id2biotype[, paste(transcript_biotype, collapse=';'), by="ensembl_gene_id_version"]
id2biotype[, gene := ensemblid2gene$GeneSymbol[match(ensembl_gene_id_version, ensemblid2gene$Geneid)]]
setnames(id2biotype, "V1", "biotype")

top.table_AD_fem$gene_name =id2biotype$gene[match(top.table_AD_fem$gene, id2biotype$ensembl_gene_id_version)]
