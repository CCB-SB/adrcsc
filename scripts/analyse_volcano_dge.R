suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(muscat))
suppressPackageStartupMessages(library(SingleCellExperiment))

setDTthreads(snakemake@threads)
snakemake@source("../scripts/helper.R")

#save.image()
#stop()
celltype_annot = snakemake@params$celltype_annot


seurat = readRDS(snakemake@input$obj)
seurat = DietSeurat(seurat, assays="RNA")
sce = as.SingleCellExperiment(seurat)
metadata = fread(snakemake@input$metadata)
pat_tbl = metadata[RVisit == 1]

sce = sce[, sce$Sample %in% pat_tbl$Sample]
sce = sce[, pat_tbl[match(sce$Sample, Sample)]$to_exclude_reason != "Healthy patient that develops disease"]

if(!celltype_annot %in% colnames(sce)){
  celltype_annot = gsub("_w_batch", "", gsub("_w_batch_and_sex", "", gsub("\\.", "-", celltype_annot)))
}


compute_genes_to_keep_per_celltype = function(sce, celltype, min_ratio_genes_expr_per_sample=0.05, min_ratio_biogroup=0.25){
  result = lapply(setdiff(unique(sce[[celltype]]), "Unknown/doublets"), function(c) {
    print(sprintf("Processing %s", c))
    sce_sub = sce[, sce[[celltype]] == c]
    gene_expressed_in_cells = pblapply(unique(sce$Sample), function(s) {
      rowMeans(counts(sce_sub[,sce_sub$Sample == s]) > 0, na.rm = T)
    })
    
    gene_expressed_in_cells_df = setDT(gene_expressed_in_cells)
    gene_expressed_in_cells_t = transpose(gene_expressed_in_cells_df)
    colnames(gene_expressed_in_cells_t) = rownames(sce)
    gene_expressed_in_cells_t[, biogroup:=name2short[as.character(sce$Diagnosis[match(unique(sce$Sample), sce$Sample)])]]
    
    gene_expressed_in_cells_per_biogroup = gene_expressed_in_cells_t[, pblapply(.SD, function(x) mean(x >= min_ratio_genes_expr_per_sample, na.rm=T)), by=biogroup]
    genes_to_consider = gene_expressed_in_cells_per_biogroup[, 2:ncol(gene_expressed_in_cells_per_biogroup)] >= min_ratio_biogroup
    rownames(genes_to_consider) = gene_expressed_in_cells_per_biogroup$biogroup

    genes_to_consider_t = as.data.frame(t(genes_to_consider))
    genes_to_consider_t[, "CIvsHC"] = genes_to_consider_t[, "HC"] | genes_to_consider_t[, "MCI"] | genes_to_consider_t[, "AD"] | genes_to_consider_t[, "PD-MCI"]
    genes_to_consider_t[, "NDvsHC"] = genes_to_consider_t[, "HC"] | genes_to_consider_t[, "MCI"] | genes_to_consider_t[, "AD"] | genes_to_consider_t[, "PD-MCI"] | genes_to_consider_t[, "PD"]
    genes_to_consider_t[, "ADvsHC"] = genes_to_consider_t[, "HC"] | genes_to_consider_t[, "AD"]
    genes_to_consider_t[, "PDvsHC"] = genes_to_consider_t[, "HC"] | genes_to_consider_t[, "PD"]
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
  names(result) = setdiff(unique(sce[[celltype]]), "Unknown/doublets")
  return(result)
}

compute_gene_per_sample_expression_stats = function(sce, celltype){
  result = lapply(setdiff(unique(sce[[celltype]]), "Unknown/doublets"), function(c) {
    print(sprintf("Processing %s", c))
    sce_sub = sce[, sce[[celltype]] == c]
    gene_mean_expressed_in_cells = pblapply(unique(sce_sub$Sample), function(s) {
      rowMeans(counts(sce_sub[,sce_sub$Sample == s]) > 0, na.rm = T)
    })
    
    aggregate_per_biogroup = function(df, agg_fun) {
      df_t = transpose(setDT(df))
      colnames(df_t) = rownames(sce_sub)
      df_t[, biogroup:=name2short[as.character(sce_sub$Diagnosis[match(unique(sce_sub$Sample), sce_sub$Sample)])]]
      df_per_biogroup = df_t[, pblapply(.SD, function(x) agg_fun(x, na.rm=T)), by=biogroup]
      df_per_biogroup_mat = as.matrix(df_per_biogroup, rownames="biogroup")
      rnames = rownames(df_per_biogroup_mat)
      df_per_biogroup_mat = apply(df_per_biogroup_mat, 2, as.numeric)
      rownames(df_per_biogroup_mat) = rnames
      return(t(df_per_biogroup_mat))
    }
    
    gene_mean_expressed_in_cells_per_biogroup_mat = aggregate_per_biogroup(gene_mean_expressed_in_cells, mean)
    gene_sd_expressed_in_cells_per_biogroup_mat = aggregate_per_biogroup(gene_mean_expressed_in_cells, sd)
    
    gene_total_expressed_in_cells = pblapply(unique(sce_sub$Sample), function(s) {
      rowSums(counts(sce_sub[,sce_sub$Sample == s]) > 0, na.rm = T)
    })
    
    gene_expressed_in_cells_per_biogroup_mat = aggregate_per_biogroup(gene_total_expressed_in_cells, sum)
    
    return(list(num_cells=gene_expressed_in_cells_per_biogroup_mat,
                mean_cells=gene_mean_expressed_in_cells_per_biogroup_mat,
                sd_cells=gene_sd_expressed_in_cells_per_biogroup_mat
                ))
  })
  names(result) = setdiff(unique(sce[[celltype]]), "Unknown/doublets")
  return(result)
}

ctbio05c005_file = sprintf("results/processed/genes_to_keep_per_%s_bio0.5_c0.05.rds", celltype_annot)
if(file.exists(ctbio05c005_file)) {
  genes_to_keep_per_celltype_bio0.5_c0.05 = readRDS(ctbio05c005_file)
} else {
  genes_to_keep_per_celltype_bio0.5_c0.05 = compute_genes_to_keep_per_celltype(sce, celltype_annot, min_ratio_genes_expr_per_sample = 0.05, min_ratio_biogroup = 0.5)
  saveRDS(genes_to_keep_per_celltype_bio0.5_c0.05, ctbio05c005_file)
}

ctbio025c005_file = sprintf("results/processed/genes_to_keep_per_%s_bio0.25_c0.05.rds", celltype_annot)
if(file.exists(ctbio025c005_file)) {
  genes_to_keep_per_celltype_bio0.25_c0.05 = readRDS(ctbio025c005_file)
} else {
  genes_to_keep_per_celltype_bio0.25_c0.05 = compute_genes_to_keep_per_celltype(sce, celltype_annot, min_ratio_genes_expr_per_sample = 0.05, min_ratio_biogroup = 0.25)
  saveRDS(genes_to_keep_per_celltype_bio0.25_c0.05, ctbio025c005_file)
}

ct_stats_file = sprintf("results/processed/gene_per_sample_stats_%s.rds", celltype_annot)
if(file.exists(ct_stats_file)) {
  gene_per_sample_stats = readRDS(ct_stats_file)
} else {
  gene_per_sample_stats = compute_gene_per_sample_expression_stats(sce, celltype_annot)
  saveRDS(gene_per_sample_stats, ct_stats_file)
}

# add columns: #AD cells, #PD cells, #MCI cells, #PD-MCI cells, #HC cells
# add columns: AD mean sample cell ratio, AD sd cell ratio, PD mean cell ratio...
# add columns: #AD samples > 0.05 ratio
save.image(celltype_annot)

dge = readRDS(snakemake@input$deg_all_col)
dge = dge[!gene %in% c("Metazoa-SRP", "Y-RNA")]

comparisons = c("ADvsHC", "PDvsHC", "MCIvsHC", "PDMCIvsHC", "CIvsHC", "NDvsHC", "ADvsPD", "AD_PDvsHC", "MCI_PDMCIvsHC", "ADvsMCI", "PDvsMCI", "AD_PDvsMCI", "ADvsMCI_PDMCI", "PDvsMCI_PDMCI", "AD_PDvsMCI_PDMCI")

pb_results = as.data.table(melt(dge[,c("gene", "cluster_id", sprintf("logFC__%s", comparisons)), with=F]))
pb_results[, contrast:=gsub("logFC__", "", variable)]
setnames(pb_results, "value", "logFC")
pb_results[, variable:=NULL]
pb_results2 = as.data.table(melt(dge[,c("gene", "cluster_id", sprintf("p_adj.loc__%s", comparisons)), with=F]))
pb_results2[, contrast:=gsub("p_adj.loc__", "", variable)]
setnames(pb_results2, "value", "p_adj.loc")
pb_results2[, variable:=NULL]

pb_results = merge(pb_results, pb_results2)

biogroups = name2short[as.character(unique(sce$Diagnosis))]

samples_per_biogroup = table(name2short[as.character(sce$Diagnosis[match(unique(sce$Sample), sce$Sample)])])

for(c in unique(pb_results$cluster_id)) {
  for(b in biogroups) {
    pb_results[cluster_id == c, (sprintf("#%s cells", b)):= gene_per_sample_stats[[c]]$num_cells[gene, b]]
    pb_results[cluster_id == c, (sprintf("%s mean sample cell ratio", b)):= gene_per_sample_stats[[c]]$mean_cells[gene, b]]
    pb_results[cluster_id == c, (sprintf("%s sd sample cell ratio", b)):= gene_per_sample_stats[[c]]$sd_cells[gene, b]]
  }
  pb_results[cluster_id == c, `#CI cells`:= gene_per_sample_stats[[c]]$num_cells[gene, "AD"] + 
                                     gene_per_sample_stats[[c]]$num_cells[gene, "MCI"] +
                                     gene_per_sample_stats[[c]]$num_cells[gene, "PD-MCI"]
      ]
  pb_results[cluster_id == c, `#ND cells`:= gene_per_sample_stats[[c]]$num_cells[gene, "AD"] + 
        gene_per_sample_stats[[c]]$num_cells[gene, "MCI"] +
        gene_per_sample_stats[[c]]$num_cells[gene, "PD-MCI"] +
        gene_per_sample_stats[[c]]$num_cells[gene, "PD"]
      ]
 for(comp in unique(pb_results$contrast)) {
   pb_results[cluster_id == c & contrast == comp, gene_expressed_keep_bio0.5_c0.05:=genes_to_keep_per_celltype_bio0.5_c0.05[[c]][gene,][[comp]]]
   pb_results[cluster_id == c & contrast == comp, gene_expressed_keep_bio0.25_c0.05:=genes_to_keep_per_celltype_bio0.25_c0.05[[c]][gene,][[comp]]]
 }
}

rmax = function(x) {
  apply(x, 1, max)
}

pb_results[contrast == "CIvsHC", mean_sample_cell_ratio:=rmax(data.frame(CI=(`AD mean sample cell ratio`*samples_per_biogroup[["AD"]] +
                                                                           `MCI mean sample cell ratio`*samples_per_biogroup[["MCI"]] +
                                                                           `PD-MCI mean sample cell ratio`*samples_per_biogroup[["PD-MCI"]]) / 
                                                                       (samples_per_biogroup[["AD"]]+samples_per_biogroup[["MCI"]]+samples_per_biogroup[["PD-MCI"]]),
                                                                     HC=`HC mean sample cell ratio`))]
pb_results[contrast == "NDvsHC", mean_sample_cell_ratio:=rmax(data.frame(ND=(`AD mean sample cell ratio`*samples_per_biogroup[["AD"]] +
                                                                            `MCI mean sample cell ratio`*samples_per_biogroup[["MCI"]] +
                                                                            `PD-MCI mean sample cell ratio`*samples_per_biogroup[["PD-MCI"]] +
                                                                            `PD mean sample cell ratio`*samples_per_biogroup[["PD"]]) / 
                                                                        (samples_per_biogroup[["AD"]]+samples_per_biogroup[["MCI"]]+samples_per_biogroup[["PD-MCI"]]+samples_per_biogroup[["PD"]]),
                                                                      HC=`HC mean sample cell ratio`))]
pb_results[contrast == "ADvsHC", mean_sample_cell_ratio:=rmax(data.frame(AD=`AD mean sample cell ratio`,
                                                                      HC=`HC mean sample cell ratio`
))]
pb_results[contrast == "PDvsHC", mean_sample_cell_ratio:=rmax(data.frame(PD=`PD mean sample cell ratio`,
                                                                      HC=`HC mean sample cell ratio`
))]
pb_results[contrast == "PDMCIvsHC", mean_sample_cell_ratio:=rmax(data.frame(PDMCI=`PD-MCI mean sample cell ratio`,
                                                                  HC=`HC mean sample cell ratio`
))]
pb_results[contrast == "MCIvsHC", mean_sample_cell_ratio:=rmax(data.frame(MCI=`MCI mean sample cell ratio`,
                                                                     HC=`HC mean sample cell ratio`
))]
pb_results[contrast == "ADvsPD", mean_sample_cell_ratio:=rmax(data.frame(`AD mean sample cell ratio`,
                                                                      `PD mean sample cell ratio`
))]

pb_results[contrast == "AD_PDvsHC", mean_sample_cell_ratio:=rmax(data.frame(AD_PD=(`AD mean sample cell ratio`*samples_per_biogroup[["AD"]] +
                                                                            `PD mean sample cell ratio`*samples_per_biogroup[["PD"]]) / (samples_per_biogroup[["AD"]] + samples_per_biogroup[["PD"]]),
                                                                  HC=`HC mean sample cell ratio`
))]
pb_results[contrast == "AD_PDvsMCI", mean_sample_cell_ratio:=rmax(data.frame(AD_PD=(`AD mean sample cell ratio`*samples_per_biogroup[["AD"]] +
                                                                            `PD mean sample cell ratio`*samples_per_biogroup[["PD"]]) / (samples_per_biogroup[["AD"]] + samples_per_biogroup[["PD"]]),
                                                                  MCI=`MCI mean sample cell ratio`
))]
pb_results[contrast == "AD_PDvsMCI_PDMCI", mean_sample_cell_ratio:=rmax(data.frame(AD_PD=(`AD mean sample cell ratio`*samples_per_biogroup[["AD"]] +
                                                                            `PD mean sample cell ratio`*samples_per_biogroup[["PD"]]) / (samples_per_biogroup[["AD"]] + samples_per_biogroup[["PD"]]),
                                                                  MCI_PDMCI=(`MCI mean sample cell ratio`*samples_per_biogroup[["MCI"]] + `PD-MCI mean sample cell ratio`*samples_per_biogroup[["PD-MCI"]]) / 
                                                                    (samples_per_biogroup[["MCI"]] + samples_per_biogroup[["PD-MCI"]])
))]
pb_results[contrast == "ADvsMCI_PDMCI", mean_sample_cell_ratio:=rmax(data.frame(AD=`AD mean sample cell ratio`,
                                                                  MCI_PDMCI=(`MCI mean sample cell ratio`*samples_per_biogroup[["MCI"]] + `PD-MCI mean sample cell ratio`*samples_per_biogroup[["PD-MCI"]]) / 
                                                                    (samples_per_biogroup[["MCI"]] + samples_per_biogroup[["PD-MCI"]])
))]
pb_results[contrast == "PDvsMCI_PDMCI", mean_sample_cell_ratio:=rmax(data.frame(PD=`PD mean sample cell ratio`,
                                                                  MCI_PDMCI=(`MCI mean sample cell ratio`*samples_per_biogroup[["MCI"]] + `PD-MCI mean sample cell ratio`*samples_per_biogroup[["PD-MCI"]]) / 
                                                                    (samples_per_biogroup[["MCI"]] + samples_per_biogroup[["PD-MCI"]])
))]

pb_results[contrast == "MCI_PDMCIvsHC", mean_sample_cell_ratio:=rmax(data.frame(MCI_PDMCI=(`MCI mean sample cell ratio`*samples_per_biogroup[["MCI"]] + `PD-MCI mean sample cell ratio`*samples_per_biogroup[["PD-MCI"]]) / 
                                                                    (samples_per_biogroup[["MCI"]] + samples_per_biogroup[["PD-MCI"]]),
                                                                  HC=`HC mean sample cell ratio`
))]
pb_results[contrast == "PDvsMCI", mean_sample_cell_ratio:=rmax(data.frame(PD=`PD mean sample cell ratio`,
                                                                  MCI=`MCI mean sample cell ratio`
))]
pb_results[contrast == "ADvsMCI", mean_sample_cell_ratio:=rmax(data.frame(AD=`AD mean sample cell ratio`,
                                                                  MCI=`MCI mean sample cell ratio`
))]

setkey(pb_results, p_adj.loc, logFC)

httr::set_config(httr::config(ssl_verifypeer=FALSE))
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

ensemblid2gene = fread(snakemake@input$geneid2symbol)
pb_results[, gene_id:=ensemblid2gene$Geneid[match(gene, ensemblid2gene$GeneSymbol)]]

id2biotype = data.table(getBM(c("ensembl_gene_id_version", "transcript_biotype"), filters="ensembl_gene_id_version", values=unique(pb_results$gene_id), mart=ensembl))
id2biotype = id2biotype[, paste(transcript_biotype, collapse=';'), by="ensembl_gene_id_version"]
id2biotype[, gene := ensemblid2gene$GeneSymbol[match(ensembl_gene_id_version, ensemblid2gene$Geneid)]]
setnames(id2biotype, "V1", "biotype")

pb_results[, biotype:=id2biotype$biotype[match(gene, id2biotype$gene)]]
pb_results[, show_label:=!grepl("pseudogene", biotype) & !grepl("^A[CFLP]\\d+\\.\\d$", gene) & !grepl("^Z\\d+\\.\\d$", gene) & !grepl("^C\\d+orf\\d+", gene) & !grepl("^RP[SL][[:digit:]]", gene)]
pb_results[, category:="Not deregulated"]
pb_results[logFC >= 0.5 & p_adj.loc <= 0.05, category:="Sig. upregulated"]
pb_results[logFC <= -0.5 & p_adj.loc <= 0.05, category:="Sig. downregulated"]

saveRDS(pb_results, snakemake@output$deg_all_col_w_add_info)
fwrite(pb_results, snakemake@output$deg_all_col_w_add_info_csv, sep='\t')
fwrite(pb_results[p_adj.loc <= 0.05], snakemake@output$deg_all_col_w_add_info_sig_p_csv, sep='\t')

df_for_labels = rbind(pb_results[(show_label) & logFC >= 0.5 & p_adj.loc <= 0.05 & mean_sample_cell_ratio >= 0.05, head(.SD, 4), by=c("contrast", "cluster_id")],
                      pb_results[(show_label) & logFC <= -0.5 & p_adj.loc <= 0.05 & mean_sample_cell_ratio >= 0.05, head(.SD, 4), by=c("contrast", "cluster_id")])

make_volcano_plot = function(dge, label_df, fc_limit=0.5, show_label=T) {
  if(!show_label){
    label_df = label_df[FALSE,]
  }
  ggplot(dge[order(abs(logFC) * -log10(p_adj.loc))], aes(x=logFC, y=-log10(p_adj.loc), fill=category, color=category, size=mean_sample_cell_ratio)) + facet_grid(rows=vars(contrast), cols = vars(cluster_id)) +
    geom_point(aes(size=mean_sample_cell_ratio), alpha=0.7) +
    scale_fill_manual(name="", values=c("Not deregulated"="lightgrey", "Sig. upregulated"="#d6604d", "Sig. downregulated"="#4393c3")) +
    scale_color_manual(name="", values=c("Not deregulated"="lightgrey", "Sig. upregulated"="#d6604d", "Sig. downregulated"="#4393c3")) +
    scale_size() + 
    geom_hline(yintercept = -log10(0.05), color="lightgrey", linetype="dashed") +
    geom_vline(xintercept = -fc_limit, color="lightgrey", linetype="dashed") + geom_vline(xintercept = fc_limit, color="lightgrey", linetype="dashed") +
    geom_text_repel(aes(x=logFC, y=-log10(p_adj.loc), label=gene), data = label_df, color="black", size=2) + 
    xlab(bquote(log[2]~fold~change)) + ylab(bquote(-log[10]~adjusted~P-value)) +
    theme_cowplot(10) + theme(legend.position = "bottom", legend.direction = "horizontal")
}

num_celltypes = length(unique(pb_results$cluster_id))
full_volcano = make_volcano_plot(pb_results, df_for_labels)
save_plot(snakemake@output$unfiltered_volcano, full_volcano, base_width=42*num_celltypes, base_height=210, unit="mm", limitsize = FALSE)
full_volcano = make_volcano_plot(pb_results, df_for_labels, show_label=F)
save_plot(snakemake@output$unfiltered_volcano_unlabeled, full_volcano, base_width=42*num_celltypes, base_height=210, unit="mm", limitsize = FALSE)

df_for_labels.stringent_bio0.25_c0.05 = rbind(pb_results[gene_expressed_keep_bio0.25_c0.05 & show_label & logFC >= 0.5 & p_adj.loc <= 0.05, head(.SD, 4), by=c("contrast", "cluster_id")],
                                              pb_results[gene_expressed_keep_bio0.25_c0.05 & show_label & logFC <= -0.5 & p_adj.loc <= 0.05, head(.SD, 4), by=c("contrast", "cluster_id")])
df_for_labels.stringent_bio0.5_c0.05 = rbind(pb_results[gene_expressed_keep_bio0.5_c0.05 & show_label & logFC >= 0.5 & p_adj.loc <= 0.05, head(.SD, 4), by=c("contrast", "cluster_id")],
                                             pb_results[gene_expressed_keep_bio0.5_c0.05 & show_label & logFC <= -0.5 & p_adj.loc <= 0.05, head(.SD, 4), by=c("contrast", "cluster_id")])

volcano_bio0.25_c0.05 = make_volcano_plot(pb_results[(gene_expressed_keep_bio0.25_c0.05)], df_for_labels.stringent_bio0.25_c0.05)
save_plot(snakemake@output$filtered_bio025_c005, volcano_bio0.25_c0.05, base_width=42*num_celltypes, base_height=210, unit="mm", limitsize = FALSE)
volcano_bio0.5_c0.05 = make_volcano_plot(pb_results[(gene_expressed_keep_bio0.5_c0.05)], df_for_labels.stringent_bio0.5_c0.05)
save_plot(snakemake@output$filtered_bio05_c005, volcano_bio0.5_c0.05, base_width=42*num_celltypes, base_height=210, unit="mm", limitsize = FALSE)

volcano_bio0.25_c0.05 = make_volcano_plot(pb_results[(gene_expressed_keep_bio0.25_c0.05)], df_for_labels.stringent_bio0.25_c0.05, show_label=F)
save_plot(snakemake@output$filtered_bio025_c005_unlabeled, volcano_bio0.25_c0.05, base_width=42*num_celltypes, base_height=210, unit="mm", limitsize = FALSE)
volcano_bio0.5_c0.05 = make_volcano_plot(pb_results[(gene_expressed_keep_bio0.5_c0.05)], df_for_labels.stringent_bio0.5_c0.05, show_label=F)
save_plot(snakemake@output$filtered_bio05_c005_unlabeled, volcano_bio0.5_c0.05, base_width=42*num_celltypes, base_height=210, unit="mm", limitsize = FALSE)
