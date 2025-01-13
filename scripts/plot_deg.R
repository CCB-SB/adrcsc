suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(superheat))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(ggplotify))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(Seurat))

setDTthreads(snakemake@threads)

snakemake@source("../scripts/shift_legend.R")
snakemake@source("../scripts/helper.R")
#save.image()
#stop()

colors = fread("data/colors.csv", strip.white = FALSE)
color_v = colors$Color
names(color_v) = colors$ID

ds_pb_edgeR.all_col = readRDS(snakemake@input$deg_all_col_w_add_info)

ds_pb_edgeR.filtered_col = ds_pb_edgeR.all_col[abs(logFC) >= 0.5 & p_adj.loc < 0.05]
sig_diff_deg = ds_pb_edgeR.filtered_col[, .N, by=c("cluster_id", "contrast")]
ds_pb_edgeR.filtered_col.stringent_bio0.25_c0.05 = ds_pb_edgeR.all_col[abs(logFC) >= 0.5 & p_adj.loc < 0.05 & gene_expressed_keep_bio0.25_c0.05]
sig_diff_deg.stringent_bio0.25_c0.05 = ds_pb_edgeR.filtered_col.stringent_bio0.25_c0.05[, .N, by=c("cluster_id", "contrast")]
ds_pb_edgeR.filtered_col.stringent_bio0.5_c0.05 = ds_pb_edgeR.all_col[abs(logFC) >= 0.5 & p_adj.loc < 0.05 & gene_expressed_keep_bio0.5_c0.05]
sig_diff_deg.stringent_bio0.5_c0.05 = ds_pb_edgeR.filtered_col.stringent_bio0.5_c0.05[, .N, by=c("cluster_id", "contrast")]

ds_pb_edgeR.filtered_col_up = ds_pb_edgeR.all_col[logFC >= 0.5 & p_adj.loc < 0.05]
sig_diff_up = ds_pb_edgeR.filtered_col_up[, .N, by=c("cluster_id", "contrast")]
ds_pb_edgeR.filtered_col_up.stringent_bio0.25_c0.05 = ds_pb_edgeR.all_col[logFC >= 0.5 & p_adj.loc < 0.05 & gene_expressed_keep_bio0.25_c0.05]
sig_diff_up.stringent_bio0.25_c0.05 = ds_pb_edgeR.filtered_col_up.stringent_bio0.25_c0.05[, .N, by=c("cluster_id", "contrast")]
ds_pb_edgeR.filtered_col_up.stringent_bio0.5_c0.05 = ds_pb_edgeR.all_col[logFC >= 0.5 & p_adj.loc < 0.05 & gene_expressed_keep_bio0.5_c0.05]
sig_diff_up.stringent_bio0.5_c0.05 = ds_pb_edgeR.filtered_col_up.stringent_bio0.25_c0.05[, .N, by=c("cluster_id", "contrast")]

ds_pb_edgeR.filtered_col_down = ds_pb_edgeR.all_col[logFC <= -0.5 & p_adj.loc < 0.05]
sig_diff_down = ds_pb_edgeR.filtered_col_down[, .N, by=c("cluster_id", "contrast")]
ds_pb_edgeR.filtered_col_down.stringent_bio0.25_c0.05 = ds_pb_edgeR.all_col[logFC <= -0.5 & p_adj.loc < 0.05 & gene_expressed_keep_bio0.25_c0.05]
sig_diff_down.stringent_bio0.25_c0.05 = ds_pb_edgeR.filtered_col_down.stringent_bio0.25_c0.05[, .N, by=c("cluster_id", "contrast")]
ds_pb_edgeR.filtered_col_down.stringent_bio0.5_c0.05 = ds_pb_edgeR.all_col[logFC <= -0.5 & p_adj.loc < 0.05 & gene_expressed_keep_bio0.5_c0.05]
sig_diff_down.stringent_bio0.5_c0.05 = ds_pb_edgeR.filtered_col_down.stringent_bio0.5_c0.05[, .N, by=c("cluster_id", "contrast")]

get_sig_diff_heatmap = function(sig_diff){
  sig_diff_heat = dcast(sig_diff, cluster_id~contrast)
  sig_diff_heat_mtx = as.matrix(sig_diff_heat[,2:ncol(sig_diff_heat)])
  rownames(sig_diff_heat_mtx) = sig_diff_heat$cluster_id
  sig_diff_heat_mtx = t(sig_diff_heat_mtx)
  comparison_order = c("NDvsHC", "CIvsHC", "ADvsHC", "PDvsHC", "PDMCIvsHC", "MCIvsHC", "ADvsPD")
  for(c in comparison_order){
    if(!c %in% rownames(sig_diff_heat_mtx)){
      sig_diff_heat_mtx = rbind(sig_diff_heat_mtx, rep(0, ncol(sig_diff_heat_mtx)))
      rownames(sig_diff_heat_mtx)[nrow(sig_diff_heat_mtx)] = c
    }
  }
  sig_diff_heat_mtx = sig_diff_heat_mtx[comparison_order,, drop=F]
  sig_diff_heat_mtx[is.na(sig_diff_heat_mtx)] = 0
  sig_diff_heat_mtx
}
  

sig_diff_heat_all = get_sig_diff_heatmap(sig_diff_deg)
sig_diff_heat_all.stringent_bio0.25_c0.05 = get_sig_diff_heatmap(sig_diff_deg.stringent_bio0.25_c0.05)
sig_diff_heat_all.stringent_bio0.5_c0.05 = get_sig_diff_heatmap(sig_diff_deg.stringent_bio0.5_c0.05)
sig_diff_heat_up = get_sig_diff_heatmap(sig_diff_up)
sig_diff_heat_up.stringent_bio0.25_c0.05 = get_sig_diff_heatmap(sig_diff_up.stringent_bio0.25_c0.05)
sig_diff_heat_up.stringent_bio0.5_c0.05 = get_sig_diff_heatmap(sig_diff_up.stringent_bio0.5_c0.05)
sig_diff_heat_down = get_sig_diff_heatmap(sig_diff_down)
sig_diff_heat_down.stringent_bio0.25_c0.05 = get_sig_diff_heatmap(sig_diff_down.stringent_bio0.25_c0.05)
sig_diff_heat_down.stringent_bio0.5_c0.05 = get_sig_diff_heatmap(sig_diff_down.stringent_bio0.5_c0.05)

scale = function (x, rows, columns) {
  if(rows){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    x = (x - m)/s
  }
  if(columns){
    x = t(x)
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    x = (x - m)/s
    x = t(x)
  }
  return(x)
}

get_text_color = function(mtcars) {
  # set the text colors 
  # identify all scaled values that fall below -0.3
  mtcars.col <- mtcars < -0.3
  # set all values that satisfy the condition to "white"
  mtcars.col <- gsub("TRUE", "white", mtcars.col)
  # set all values that do not satisfy the condition to "black"
  mtcars.col <- gsub("FALSE", "black", mtcars.col)
  # convert to matrix
  mtcars.col <- matrix(mtcars.col, ncol = ncol(mtcars))
  return(mtcars.col)
}

plot_deg_heatmap = function(heat, filename){
  if(ncol(heat) == 1){
    pheatmap(sig_diff_heat_all, color = scales::gradient_n_pal(rev(brewer.pal(9, "RdYlBu")))(seq(0,1,length.out=100)),
             cluster_rows=F, cluster_cols=F, display_numbers = T,
             number_format = "%d", number_color = "black",
             angle_col = 0, border_color = "white", legend = F,
             show_rownames = T, show_colnames = T,
             filename = filename, width=5.51181, height=3.34646, unit="mm")
  } else {
    plot = superheat(heat, col.dendrogram = F, row.dendrogram = F,
              X.text=heat, scale=F,
              grid.hline.col = "white", grid.vline.col = "white",
              legend = F,
              legend.text.size = 8,
              X.text.size = 3,
              X.text.col = get_text_color(heat),
              bottom.label.text.size = 3, left.label.text.size = 3,
              bottom.label.text.angle = 90,
              heat.pal = rev(brewer.pal(9, "RdYlBu"))
    )
    save_plot(filename, plot$plot, base_width=140, base_height = 85, unit = "mm")
  }
}

plot_deg_heatmap(sig_diff_heat_all, snakemake@output$all_comp_heat)
plot_deg_heatmap(sig_diff_heat_all.stringent_bio0.25_c0.05, snakemake@output$all_comp_heat_bio025_c005)
plot_deg_heatmap(sig_diff_heat_all.stringent_bio0.5_c0.05, snakemake@output$all_comp_heat_bio05_c005)
plot_deg_heatmap(sig_diff_heat_up, snakemake@output$up_comp_heat)
plot_deg_heatmap(sig_diff_heat_up.stringent_bio0.25_c0.05, snakemake@output$up_comp_heat_bio025_c005)
plot_deg_heatmap(sig_diff_heat_up.stringent_bio0.5_c0.05, snakemake@output$up_comp_heat_bio05_c005)
plot_deg_heatmap(sig_diff_heat_down, snakemake@output$down_comp_heat)
plot_deg_heatmap(sig_diff_heat_down.stringent_bio0.25_c0.05, snakemake@output$down_comp_heat_bio025_c005)
plot_deg_heatmap(sig_diff_heat_down.stringent_bio0.5_c0.05, snakemake@output$down_comp_heat_bio05_c005)


############# upset plots #####################
############# all diseases and celltypes together
get_celltype_label_colors = function(upset_mat) {
  upset_mat = upset_mat[,colSums(upset_mat) > 0, drop=F]
  comp_label_order = colnames(upset_mat)[order(-colSums(upset_mat), colnames(upset_mat))]
  
  comp_label_ct = unlist(lapply(strsplit(comp_label_order, "__"), function(x) x[1]))
  comp_label_colors = color_v[comp_label_ct]
  comp_label_colors
}

create_upset_plot = function(df, df_up, df_down, filename_all, filename_up, filename_down) {
  get_upset_mat = function(df) {
    de_genes_per_celltype_and_contrast = list()
    for(c in unique(df$cluster_id)) {
      for(comp in unique(df$contrast)) {
        de_genes_per_celltype_and_contrast[[sprintf("%s__%s", c, comp)]] = df[cluster_id == c & contrast == comp]$gene
      }
    }
    fromList(de_genes_per_celltype_and_contrast)
  }
  upset_mat = get_upset_mat(df)
  upset_mat = upset_mat[, colSums(upset_mat) > 0, drop=F] # drop comparisons that have no genes
  comp_label_colors = get_celltype_label_colors(upset_mat)
  comp_label_colors = comp_label_colors[1:min(30, length(comp_label_colors))]
  
  upset_mat_up = get_upset_mat(df_up)
  for(c in colnames(upset_mat)[!colnames(upset_mat) %in% colnames(upset_mat_up)]){
    upset_mat_up[[c]] = 0
  }
  upset_mat_up = upset_mat_up[,colnames(upset_mat)]
  
  upset_mat_down = get_upset_mat(df_down)
  for(c in colnames(upset_mat)[!colnames(upset_mat) %in% colnames(upset_mat_down)]){
    upset_mat_down[[c]] = 0
  }
  upset_mat_down = upset_mat_down[,colnames(upset_mat)]

  sets = colnames(upset_mat)[order(-colSums(upset_mat), colnames(upset_mat))]
  sets = sets[1:min(30, length(sets))]
  
  tryCatch({
    p = UpSetR::upset(upset_mat, order.by="freq", keep.order = T, sets=sets, nsets=length(comp_label_colors), nintersects = 40, mb.ratio=c(0.5,0.5), sets.bar.color = comp_label_colors)
    layout = c(area(t=1, l=2, b=2, r=31),
               area(t=2, l=1, b=2, r=7))
    plot = (as.ggplot(p) + (as.ggplot(p$Sizes))) + plot_layout(design=layout)
    save_plot(filename_all, plot, base_width=270, base_height = 200, unit = "mm")
  },
  error=function(x){
    file.create(filename_all)
  })
  
  tryCatch({
    p = UpSetR::upset(upset_mat_up, order.by="freq", keep.order = T, sets=sets, nsets=length(comp_label_colors), nintersects = 40, mb.ratio=c(0.5,0.5), sets.bar.color = comp_label_colors)
    layout = c(area(t=1, l=2, b=2, r=31),
               area(t=2, l=1, b=2, r=7))
    plot = (as.ggplot(p) + (as.ggplot(p$Sizes))) + plot_layout(design=layout)
    save_plot(filename_up, plot, base_width=270, base_height = 200, unit = "mm")
  },
  error=function(x){
    file.create(filename_up)
  })
  
  tryCatch({
    p = UpSetR::upset(upset_mat_down, order.by="freq", keep.order = T, sets=sets, nsets=length(comp_label_colors), nintersects = 40, mb.ratio=c(0.5,0.5), sets.bar.color = comp_label_colors)
    layout = c(area(t=1, l=2, b=2, r=31),
               area(t=2, l=1, b=2, r=7))
    plot = (as.ggplot(p) + (as.ggplot(p$Sizes))) + plot_layout(design=layout)
    save_plot(filename_down, plot, base_width=270, base_height = 200, unit = "mm")
  },
  error=function(x){
    file.create(filename_down)
  })
}

create_upset_plot(ds_pb_edgeR.filtered_col, ds_pb_edgeR.filtered_col_up, ds_pb_edgeR.filtered_col_down, snakemake@output$all_comp_upset, snakemake@output$up_comp_upset, snakemake@output$down_comp_upset)
create_upset_plot(ds_pb_edgeR.filtered_col.stringent_bio0.25_c0.05, ds_pb_edgeR.filtered_col_up.stringent_bio0.25_c0.05, ds_pb_edgeR.filtered_col_down.stringent_bio0.25_c0.05, snakemake@output$all_comp_upset_bio025_c005, snakemake@output$up_comp_upset_bio025_c005, snakemake@output$down_comp_upset_bio025_c005)
create_upset_plot(ds_pb_edgeR.filtered_col.stringent_bio0.5_c0.05, ds_pb_edgeR.filtered_col_up.stringent_bio0.5_c0.05, ds_pb_edgeR.filtered_col_down.stringent_bio0.5_c0.05, snakemake@output$all_comp_upset_bio05_c005, snakemake@output$up_comp_upset_bio05_c005, snakemake@output$down_comp_upset_bio05_c005)

############### only per disease
get_disease_comp_label_colors = function(upset_mat) {
  comp_label_order = names(sort(-colSums(upset_mat)))
  color_mapping=c("PDvsHC"=unname(color_v["PD"]),
                  "CIvsHC"=unname(color_v["Cognitive Impairment"]),
                  "NDvsHC"=unname(color_v["Neurodegeneration"]),
                  "ADvsPD"="#8c510a",
                  "ADvsHC"=unname(color_v["AD"]),
                  "MCIvsHC"=unname(color_v["MCI"]),
                  "PDMCIvsHC"=unname(color_v["PD-MCI"]))
  comp_label_colors = color_mapping[comp_label_order]
  comp_label_colors
}

plot_upset_per_disease = function(df, filename) {
  de_genes_per_celltype_and_contrast = list()
  for(c in unique(df$cluster_id)) {
    for(comp in unique(df$contrast)) {
      de_genes_per_celltype_and_contrast[[comp]] = unique(df[contrast == comp]$gene)
    }
  }
  
  fromList = function (input) 
  {
    elements <- unique(unlist(input))
    data <- unlist(lapply(input, function(x) {
      x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input), byrow = F))
    data <- data[which(rowSums(data) != 0), , drop=FALSE]
    names(data) <- names(input)
    return(data)
  }
  
  
  upset_mat = fromList(de_genes_per_celltype_and_contrast)
  upset_mat = upset_mat[, colSums(upset_mat) > 0, drop=F] # drop comparisons that have no genes
  comp_label_colors = get_disease_comp_label_colors(upset_mat)

  tryCatch({
    p = UpSetR::upset(upset_mat, order.by="freq", nsets=25, nintersects = 40, mb.ratio=c(0.5,0.5), sets.bar.color = comp_label_colors)
    layout = c(area(t=1, l=2, b=2, r=31),
               area(t=2, l=1, b=2, r=7))
    plot = (as.ggplot(p) + (as.ggplot(p$Sizes))) + plot_layout(design=layout)
    save_plot(filename, plot, base_width=270, base_height = 200, unit = "mm")
  }, error=function(x){
    file.create(filename)
  })
}

plot_upset_per_disease(ds_pb_edgeR.filtered_col, snakemake@output$disease_comp_upset)
plot_upset_per_disease(ds_pb_edgeR.filtered_col.stringent_bio0.25_c0.05, snakemake@output$disease_comp_bio025_c005_upset)
plot_upset_per_disease(ds_pb_edgeR.filtered_col.stringent_bio0.5_c0.05, snakemake@output$disease_comp_bio05_c005_upset)
plot_upset_per_disease(ds_pb_edgeR.filtered_col_down, snakemake@output$disease_down_comp_upset)
plot_upset_per_disease(ds_pb_edgeR.filtered_col_down.stringent_bio0.25_c0.05, snakemake@output$disease_down_comp_bio025_c005_upset)
plot_upset_per_disease(ds_pb_edgeR.filtered_col_down.stringent_bio0.5_c0.05, snakemake@output$disease_down_comp_bio05_c005_upset)
plot_upset_per_disease(ds_pb_edgeR.filtered_col_up, snakemake@output$disease_up_comp_upset)
plot_upset_per_disease(ds_pb_edgeR.filtered_col_up.stringent_bio0.25_c0.05, snakemake@output$disease_up_comp_bio025_c005_upset)
plot_upset_per_disease(ds_pb_edgeR.filtered_col_up.stringent_bio0.5_c0.05, snakemake@output$disease_up_comp_bio05_c005_upset)


############### only per celltype

plot_upset_per_celltype = function(df, filename){
  de_genes_per_celltype_and_contrast = list()
  for(c in unique(df$cluster_id)) {
    de_genes_per_celltype_and_contrast[[c]] = unique(df[cluster_id == c]$gene)
  }
  
  if(length(unique(df$cluster_id)) == 1){
    file.create(filename)
    return()
  }
  
  
  upset_mat = fromList(de_genes_per_celltype_and_contrast)
  upset_mat = upset_mat[, colSums(upset_mat) > 0, drop=F] # drop comparisons that have no genes
  comp_label_colors = get_celltype_label_colors(upset_mat)
  print(comp_label_colors)
  p = UpSetR::upset(upset_mat, order.by="freq", nsets=30, nintersects = 40, mb.ratio=c(0.5,0.5))#, sets.bar.color = comp_label_colors)
  layout = c(area(t=1, l=2, b=2, r=31),
             area(t=2, l=1, b=2, r=7))
  plot = (as.ggplot(p) + (as.ggplot(p$Sizes))) + plot_layout(design=layout)
  save_plot(filename, plot, base_width=270, base_height = 200, unit = "mm")
}

plot_upset_per_celltype(ds_pb_edgeR.filtered_col, snakemake@output$celltype_comp_upset)
plot_upset_per_celltype(ds_pb_edgeR.filtered_col.stringent_bio0.25_c0.05, snakemake@output$celltype_comp_bio025_c005_upset)
plot_upset_per_celltype(ds_pb_edgeR.filtered_col.stringent_bio0.5_c0.05, snakemake@output$celltype_comp_bio05_c005_upset)
plot_upset_per_celltype(ds_pb_edgeR.filtered_col_down, snakemake@output$celltype_down_comp_upset)
plot_upset_per_celltype(ds_pb_edgeR.filtered_col_down.stringent_bio0.25_c0.05, snakemake@output$celltype_down_comp_bio025_c005_upset)
plot_upset_per_celltype(ds_pb_edgeR.filtered_col_down.stringent_bio0.5_c0.05, snakemake@output$celltype_down_comp_bio05_c005_upset)
plot_upset_per_celltype(ds_pb_edgeR.filtered_col_up, snakemake@output$celltype_up_comp_upset)
plot_upset_per_celltype(ds_pb_edgeR.filtered_col_up.stringent_bio0.25_c0.05, snakemake@output$celltype_up_comp_bio025_c005_upset)
plot_upset_per_celltype(ds_pb_edgeR.filtered_col_up.stringent_bio0.5_c0.05, snakemake@output$celltype_up_comp_bio05_c005_upset)


###################

sce = as.SingleCellExperiment(readRDS(snakemake@input$obj))
celltype_annot = gsub("_w_batch", "", gsub("_w_batch_and_sex", "", snakemake@params$celltype_annot))

if(!celltype_annot %in% colnames(sce)){
  celltype_annot = gsub("\\.", "-", celltype_annot)
}

corr_eqn = function(x,y, method="pearson", digits = 2) {
  corr_coef = round(cor(x, y, method=method), digits = digits)
  eq = substitute(~~italic(r)~"="~r2, list(r2=corr_coef))
  as.character(as.expression(eq))
}

plot_correlation_deg_cells = function(sig_heat, sce, filename) {
  sig_diff_plot_df = reshape2::melt(sig_heat, value.name="Count", varnames=c("Comp", "Celltype"))
  n_celltypes = as.vector(table(sce[[celltype_annot]]))
  names(n_celltypes) = names(table(sce[[celltype_annot]]))
  sig_diff_plot_df$total_cells = n_celltypes[sig_diff_plot_df$Celltype]
  
  correlation_labels = data.frame(Comp=unique(sig_diff_plot_df$Comp))
  correlation_labels$x = unlist(lapply(correlation_labels$Comp, function(c) (max(sig_diff_plot_df$Count[sig_diff_plot_df$Comp == c]) + min(sig_diff_plot_df$Count[sig_diff_plot_df$Comp == c])) / 2))
  correlation_labels$y = 0.99 * max(sig_diff_plot_df$total_cells)
  correlation_labels$label = unlist(lapply(correlation_labels$Comp, function(c) {
    corr_eqn(sig_diff_plot_df[sig_diff_plot_df$Comp == c,]$Count, sig_diff_plot_df[sig_diff_plot_df$Comp == c,]$total_cells, "pearson")
  }))
  
  if(max(sig_diff_plot_df$total_cells) > 1e4){
    sig_diff_plot_df$total_cells = sig_diff_plot_df$total_cells / 1e3
    correlation_labels$y = correlation_labels$y / 1e3
    yl = bquote("#Cells ["~10^3~"]")
  } else {
    yl = "#Cells"
  }
  
  if(length(unique(sig_diff_plot_df$Celltype)) <= 2){
    correlation_labels = correlation_labels[F,]
  }
  
  
  p = ggplot(sig_diff_plot_df, aes(x=Count, y=total_cells)) + facet_wrap(~factor(Comp, levels=c("NDvsHC", "CIvsHC", "ADvsHC", "MCIvsHC", "PDMCIvsHC", "PDvsHC", "ADvsPD")), scales="free_x") +
    geom_point(aes(fill=Celltype, color=Celltype), size=4) +
    geom_smooth(method="lm", se = F, color="black") +
    geom_text(aes(label = label, x=x, y=y), parse=T, color="black", data=correlation_labels) +
    scale_fill_manual(values=color_v) +
    scale_color_manual(values=color_v) +
    scale_y_continuous(expand=expand_scale(mult=.1)) +
    xlab("#DEG") + ylab(yl) +
    theme_cowplot(10) + theme(legend.position = "bottom", legend.direction = "horizontal")
  
  save_plot(filename, p, base_width=200, base_height=150, unit="mm")
}

plot_correlation_deg_cells(sig_diff_heat_all, sce, snakemake@output$deg_vs_cells)
plot_correlation_deg_cells(sig_diff_heat_all.stringent_bio0.25_c0.05, sce, snakemake@output$deg_vs_cells_bio025_c005)
plot_correlation_deg_cells(sig_diff_heat_all.stringent_bio0.5_c0.05, sce, snakemake@output$deg_vs_cells_bio05_c005)


plot_deregulation_dir = function(df_up, df_down, filename_total, filename_perc) {
  sig_diff_overall_plot_df = rbind(data.table(reshape2::melt(df_down, varnames=c("Comp", "Celltype"), value="Count"), Direction="Down"), 
                                   data.table(reshape2::melt(df_up, varnames=c("Comp", "Celltype"), value="Count"), Direction="Up"))
  
  p = ggplot(sig_diff_overall_plot_df, aes(x=Comp, y=value, fill=Direction)) + 
    geom_col(position="identity") + facet_wrap(~Celltype) +
    scale_fill_npg() +
    ylab("#DEG") +
    xlab("") +
    theme_cowplot(10) + theme(axis.text.x = element_text(angle=90))
  
  n_celltypes = length(unique(sig_diff_overall_plot_df$Celltype))
  if(n_celltypes %% 4 != 0 && n_celltypes %% 3 != 0){
    p = shift_legend(p)
  }
  save_plot(filename_total, p, base_width=140, base_height = 100, unit="mm")
  
  p = ggplot(sig_diff_overall_plot_df, aes(x=Comp, y=value, fill=Direction)) + 
    geom_col(position="fill") + facet_wrap(~Celltype) +
    scale_fill_npg() +
    scale_y_continuous(labels=scales::percent_format()) +
    ylab("DEG (%)") +
    xlab("") +
    theme_cowplot(10) + theme(axis.text.x = element_text(angle=90))
  
  if(n_celltypes %% 4 != 0 && n_celltypes %% 3 != 0){
    p = shift_legend(p)
  }
  save_plot(filename_perc, p, base_width=140, base_height = 100, unit="mm")
}

plot_deregulation_dir(sig_diff_heat_up, sig_diff_heat_down, snakemake@output$deregulation_dir, snakemake@output$deregulation_dir_perc)
plot_deregulation_dir(sig_diff_heat_up.stringent_bio0.25_c0.05, sig_diff_heat_down.stringent_bio0.25_c0.05, snakemake@output$deregulation_dir_bio025_c005, snakemake@output$deregulation_dir_bio025_c005_perc)
plot_deregulation_dir(sig_diff_heat_up.stringent_bio0.5_c0.05, sig_diff_heat_down.stringent_bio0.5_c0.05, snakemake@output$deregulation_dir_bio05_c005, snakemake@output$deregulation_dir_bio05_c005_perc)
