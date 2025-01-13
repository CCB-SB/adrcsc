# Code used for computing scRNAseq Brain sample DEGs

## preprocessing.py
- merging and collecting the gathered studies into single h5ad objects
- includes scrublet + QC metrics
- write h5ad file containing the raw and processed data

## scVI_blackbox.py
- integrate scAtlas data using scVI
- intergated representation was mainly used for checking cell type labelling and clustering across studies
- the integrated representation was not used to perform any DEG computations

## detect_marker_genes.R
- computes DEGs between control and condition on subset of single cell types
- input data is split by sex into female (F) and male (M) samples

