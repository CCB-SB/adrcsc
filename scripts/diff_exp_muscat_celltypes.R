library(SingleCellExperiment)
library(Seurat)
library(muscat)
library(limma)
library(data.table)
library(pbapply)

source("scripts/helper.R")

metadata = fread(snakemake@input$metadata)
pat_tbl = metadata[RVisit == 1]

pbmc = readRDS(snakemake@input$obj)
#pbmc$celltype_cluster <- pbmc$celltype_clusters
saveRDS(pbmc, snakemake@input$obj)

set.seed(42)
pbmc = readRDS(snakemake@input$obj)
pbmc$RVisit = metadata$RVisit[match(pbmc$Sample, metadata$Sample)]
pbmc$biogroup_short = factor(make.names(name2short[as.character(pbmc$Diagnosis)]))

pbmc = DietSeurat(pbmc, assays="RNA")
pbmc = pbmc[,pbmc$celltype != "Unknown/doublets"]

pbmc_sce = as.SingleCellExperiment(pbmc)
pbmc_sce = pbmc_sce[, pbmc_sce$Sample %in% pat_tbl$Sample]
pbmc_sce = pbmc_sce[, pat_tbl[match(pbmc_sce$Sample, Sample)]$to_exclude_reason != "Healthy patient that develops disease"]
 
sce <- prepSCE(pbmc_sce,
              kid = snakemake@params$celltype,
              gid = "biogroup_short",
              sid = "Sample",
              drop = TRUE)
saveRDS(sce, snakemake@output$sce, compress=FALSE)
pb = aggregateData(sce, assay="counts", fun="sum")
saveRDS(pb, snakemake@output$pb, compress=FALSE)

pb_sub = pb[,colnames(pb)[pat_tbl[match(colnames(pb), Sample)]$to_exclude_reason != "Healthy patient that develops disease"]]

ei <- metadata(pb_sub)$experiment_info
design = model.matrix(~ 0 + ei$group_id)
dimnames(design) <- list(ei$sample_id, levels(ei$group_id))

contrast = makeContrasts(NDvsHC=(AD+MCI+PD+PD.MCI)/4-HC,
                         CIvsHC=(AD+MCI+PD.MCI)/3-HC,
                         ADvsHC=AD-HC,
                         PDvsHC=PD-HC,
                         MCIvsHC=MCI-HC,
                         PDMCIvsHC=PD.MCI-HC,
                         ADvsPD=AD-PD,
                         AD_PDvsHC=(AD+PD)/2-HC,
                         MCI_PDMCIvsHC=(MCI+PD.MCI)/2-HC,
                         ADvsMCI=AD-MCI,
                         PDvsMCI=PD-MCI,
                         AD_PDvsMCI=(AD+PD)/2-MCI,
                         ADvsMCI_PDMCI=AD-(MCI+PD.MCI)/2,
                         PDvsMCI_PDMCI=PD-(MCI+PD.MCI)/2,
                         AD_PDvsMCI_PDMCI=(AD+PD)/2-(MCI+PD.MCI)/2,
                         levels=design)

# ds_pb_edgeR = pbDS(pb_sub, method="edgeR", design=design, contrast = contrast)
# saveRDS(ds_pb_edgeR, file.path(snakemake@params$prefix, sprintf("ds_pb_edgeR_%s.rds", snakemake@params$celltype)), compress=FALSE)
ds_pb_limma_voom = pbDS(pb_sub, method="limma-voom", design=design, contrast = contrast)
saveRDS(ds_pb_limma_voom, file.path(snakemake@params$prefix, sprintf("ds_pb_limma_voom_%s.rds", snakemake@params$celltype)), compress=FALSE)

# ds_pb_deseq2 = pbDS(pb_sub, method="DESeq2", design=design, contrast = contrast)
# saveRDS(ds_pb_deseq2, file.path(snakemake@params$prefix, sprintf("ds_pb_deseq2_%s.rds", snakemake@params$celltype)), compress=FALSE)
# 
# ds_pb_limmat = pbDS(pb_sub, method="limma-trend", design=design, contrast = contrast)
# saveRDS(ds_pb_limmat, file.path(snakemake@params$prefix, sprintf("ds_pb_limma_trend_%s.rds", snakemake@params$celltype)), compress=FALSE)

genes_to_keep_per_celltype = lapply(unique(sce$cluster_id), function(c) {
  sce_sub = sce[, sce$cluster_id == c]
  gene_expressed_in_cells = pblapply(unique(sce$sample_id), function(s) {
    rowMeans(counts(sce_sub[,sce_sub$sample_id == s]) > 0, na.rm = T)
  })
  
  gene_expressed_in_cells_df = setDT(gene_expressed_in_cells)
  gene_expressed_in_cells_t = transpose(gene_expressed_in_cells_df)
  colnames(gene_expressed_in_cells_t) = rownames(sce)
  gene_expressed_in_cells_t[, biogroup:=metadata$Diagnosis[match(unique(sce$sample_id), metadata$Sample)]]
  min_expr_cells_per_genes_per_sample = 0.05
  gene_expressed_in_cells_per_biogroup = gene_expressed_in_cells_t[, pblapply(.SD, function(x) mean(x >= min_expr_cells_per_genes_per_sample, na.rm=T)), by=biogroup]
  gene_expressed_in_at_least_ratio_biogroup = 0.5
  genes_to_consider = gene_expressed_in_cells_per_biogroup[, 2:ncol(gene_expressed_in_cells_per_biogroup)] >= gene_expressed_in_at_least_ratio_biogroup
  rownames(genes_to_consider) = name2short[gene_expressed_in_cells_per_biogroup$biogroup]
  #genes_to_keep = colSums(gene_expressed_in_cells_per_biogroup[, 2:ncol(gene_expressed_in_cells_per_biogroup)] >= gene_expressed_in_at_least_ratio_biogroup) >= 1
  genes_to_consider_t = as.data.frame(t(genes_to_consider))
  genes_to_consider_t[, "CIvsHC"] = genes_to_consider_t[, "HC"] | genes_to_consider_t[, "MCI"] | genes_to_consider_t[, "AD"] | genes_to_consider_t[, "PD-MCI"]
  genes_to_consider_t[, "NDvsHC"] = genes_to_consider_t[, "HC"] | genes_to_consider_t[, "MCI"] | genes_to_consider_t[, "AD"] | genes_to_consider_t[, "PD-MCI"] | genes_to_consider_t[, "PD"]
  genes_to_consider_t[, "ADvsHC"] = genes_to_consider_t[, "HC"] | genes_to_consider_t[, "AD"]
  genes_to_consider_t[, "PDvsHC"] = genes_to_consider_t[, "HC"] | genes_to_consider_t[, "PD"] | genes_to_consider_t[, "PD-MCI"]
  genes_to_consider_t[, "ADvsPD"] = genes_to_consider_t[, "PD"] | genes_to_consider_t[, "AD"]
  genes_to_consider_t[, "MCIvsHC"] = genes_to_consider_t[, "HC"] | genes_to_consider_t[, "MCI"]
  genes_to_consider_t[, "PDMCIvsHC"] = genes_to_consider_t[, "HC"] | genes_to_consider_t[, "PD-MCI"]
  genes_to_consider_t[, "ADvsPD"] = genes_to_consider_t[, "AD"] | genes_to_consider_t[, "PD"]
  genes_to_consider_t[, "AD_PDvsHC"] = genes_to_consider_t[, "AD"] | genes_to_consider_t[, "PD"] | genes_to_consider_t[, "HC"]
  genes_to_consider_t[, "MCI_PDMCIvsHC"] = genes_to_consider_t[, "MCI"] | genes_to_consider_t[, "PD-MCI"] | genes_to_consider_t[, "HC"]
  genes_to_consider_t[, "ADvsMCI"] = genes_to_consider_t[, "AD"] | genes_to_consider_t[, "MCI"]
  genes_to_consider_t[, "PDvsMCI"] = genes_to_consider_t[, "PD"] | genes_to_consider_t[, "MCI"]
  genes_to_consider_t[, "AD_PDvsMCI"] = genes_to_consider_t[, "AD"] | genes_to_consider_t[, "PD"] | genes_to_consider_t[, "MCI"]
  genes_to_consider_t[, "ADvsMCI_PDMCI"] = genes_to_consider_t[, "AD"] | genes_to_consider_t[, "PD-MCI"] | genes_to_consider_t[, "MCI"]
  genes_to_consider_t[, "PDvsMCI_PDMCI"] = genes_to_consider_t[, "PD"] | genes_to_consider_t[, "PD-MCI"] | genes_to_consider_t[, "MCI"]
  genes_to_consider_t[, "AD_PDvsMCI_PDMCI"] = genes_to_consider_t[, "AD"] | genes_to_consider_t[, "PD"] | genes_to_consider_t[, "PD-MCI"] | genes_to_consider_t[, "MCI"]
  return(genes_to_consider_t)
})

names(genes_to_keep_per_celltype) = unique(sce$cluster_id)
saveRDS(genes_to_keep_per_celltype, snakemake@output$genes_to_keep, compress=FALSE)
# 
# ds_pb_edgeR.all_col = as.data.table(resDS(sce, ds_pb_edgeR, bind="col"))
# saveRDS(ds_pb_edgeR.all_col, file.path(snakemake@params$prefix, sprintf("ds_pb_edgeR_all_col_%s.rds", snakemake@params$celltype)), compress=FALSE)
# ds_pb_edgeR.all_row = as.data.table(resDS(sce, ds_pb_edgeR, bind="row"))
# saveRDS(ds_pb_edgeR.all_row, file.path(snakemake@params$prefix, sprintf("ds_pb_edgeR_all_row_%s.rds", snakemake@params$celltype)), compress=FALSE)
# 
ds_pb_limma_voom.all_col = as.data.table(resDS(sce, ds_pb_limma_voom, bind="col"))
saveRDS(ds_pb_limma_voom.all_col, file.path(snakemake@params$prefix, sprintf("ds_pb_limma_voom_all_col_%s.rds", snakemake@params$celltype)), compress=FALSE)
ds_pb_limma_voom.all_row = as.data.table(resDS(sce, ds_pb_limma_voom, bind="row"))
saveRDS(ds_pb_limma_voom.all_row, file.path(snakemake@params$prefix, sprintf("ds_pb_limma_voom_all_row_%s.rds", snakemake@params$celltype)), compress=FALSE)
# 
# ds_pb_limmat.all_col = as.data.table(resDS(sce, ds_pb_limmat, bind="col"))
# saveRDS(ds_pb_limmat.all_col, file.path(snakemake@params$prefix, sprintf("ds_pb_limma_trend_all_col_%s.rds", snakemake@params$celltype)), compress=FALSE)
# ds_pb_limmat.all_row = as.data.table(resDS(sce, ds_pb_limmat, bind="row"))
# saveRDS(ds_pb_limmat.all_row, file.path(snakemake@params$prefix, sprintf("ds_pb_limma_trend_all_row_%s.rds", snakemake@params$celltype)), compress=FALSE)
# 
# ds_pb_deseq2.all_col = as.data.table(resDS(sce, ds_pb_deseq2, bind="col"))
# saveRDS(ds_pb_deseq2.all_col, file.path(snakemake@params$prefix, sprintf("ds_pb_deseq2_all_col_%s.rds", snakemake@params$celltype)), compress=FALSE)
# ds_pb_deseq2.all_row = as.data.table(resDS(sce, ds_pb_deseq2, bind="row"))
# saveRDS(ds_pb_deseq2.all_row, file.path(snakemake@params$prefix, sprintf("ds_pb_deseq2_all_row_%s.rds", snakemake@params$celltype)), compress=FALSE)
