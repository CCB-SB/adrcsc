suppressPackageStartupMessages({
  library(MAST)
  library(Seurat)
  library(SeuratDisk)
  library(hdf5r)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(stringr)
  library(rsvd)
  library(scater)
  library(future)
})




compute_marker = function(seurat,id,cell_t,cond,genes,latent_var){
   """ Compute marker

    Args:
        seurat (seurat object): Seurat object  
        cell_t (string)     : Cell type to compute DEGs between condition and control
        cond   (string)     : Condition of interest vs control
        genes  (list)       : List of genes which should be tested for differential expression
        latent_var (vector) : vector of string containing the latent variable keys
    """   
  Idents(object = seurat)
  
  # subsetting seurat object to only contain cell types of interest
  sset = seurat[, seurat@meta.data[, "sub_cell_type"] == cell_t]
  sset@meta.data = seurat@meta.data[colnames(sset),]      
  
  # Split dataset by sex
  sset_m = sset[, sset@meta.data[, "sex"] == "M"]
  sset_f = sset[, sset@meta.data[, "sex"] == "F"]  
  
  sset_m <- SetIdent(sset_m, value = sset_m@meta.data$condition)
  sset_f <- SetIdent(sset_f, value = sset_f@meta.data$condition)
  
  print(paste("processing",cond,cell_t))
  
  # compute and save markers using MAST from withing seurat
  mast_markers <-FindMarkers(object =sset_m,ident.1=cond,ident.2 ="CT",max.cells.per.ident = sum(Idents(sset_m) == cond),
                           test.use = "MAST",latent.vars =latent_var,logfc.threshold = 0,min.pct = 0,features = genes)
  
  write.csv(mast_markers,paste0("./",id,"_",cond,"_ct_",cell_t,"_M_",latent_var,".csv"), row.names = TRUE)
  
  mast_markers <-FindMarkers(object =sset_f,ident.1=cond,ident.2 ="CT",max.cells.per.ident = sum(Idents(sset_f) == cond),
                             test.use = "MAST",latent.vars = latent_var,logfc.threshold = 0,min.pct = 0,features = genes)
  
  write.csv(mast_markers,paste0("./",id,"_",cond,"_ct_",cell_t,"_F_",latent_var,".csv"), row.names = TRUE)

}

plan("multiprocess", workers = 124)
# cell types which are supposed to get tested for DEG
cell_types = list("Microglia","Astrocyte",'Excitatory neuron','Inhibitory neuron',"Oligodendrocyte","OPC","Vascular","Pvalb neuron")

# loading expression data and genes to test
seurat = readRDS("/local/scRNA_brain_atlas/v1.0/h_cortex_seurat/seurat.rds")
row.names(seurat@meta.data) = seurat@meta.data$X
genes = readLines("/home/mflotho/single-cell-rna-seq-reference-atlas/ADRC_genes/all")
# defining the study-id (dataset) as latent var
latent_var = c("dataset")

for (j in 1:length(cell_types)){
    compute_marker(seurat,files[[1]],cell_types[[j]],cond = "AD",genes = genes,latent_var = latent_var)
}

packageVersion("pvca")
