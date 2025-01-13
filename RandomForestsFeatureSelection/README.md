# Feature Selection and Classification Workflow

This directory contains the workflow described in Supplementary Note 1: Machine Learning of Potential Biomarker Panels Using Random Forests. It is a Snakemake pipeline designed for feature selection in scRNA-seq data, utilizing Random Forests with backward feature elimination, as well as two deep learning-based methods (DeepLIFT and FeatureAblation) and a regression-based method (ElasticNet). The selected gene sets are evaluated using a kNN classifier and two test sets, comprising previously unseen data from independent and dependent donors. Comparing the performance on the two test sets demonstrates the critical importance of patient-level data splitting in scRNA-seq analyses to mitigate overfitting to donor-specific effects. For more information please view our paper _"A single-cell atlas to map sex-specific gene-expression changes in blood upon neurodegeneration" (Grandke et al., 2025)_.

## Prerequisites
In order to run the workflow a working Conda and Snakemake installation are required.  Additionally, scDeepFeatures is required. You can clone the official repository (https://github.com/PYangLab/scDeepFeatures) and edit the path to the the main script of scDeepFeatures (scDeepFeatures/Utils/Feature_selection_methods/Mlp/main.py) in the config.yaml file of the pipeline. The pipeline should install all other dependencies on first execution using Snakemake's Conda integration. 

## Usage
Prior to execution the config.yaml of the pipeline needs to be edited. In order to reproduce our results (or to run the workflow on your own data) the file paths for a Seurat object with training donors and a Seurat object with independent test donors needs to be set. Additionally, a csv file containing differentially expressed genes (DEGs) is required for comparison with the machine learning-based feature selections. The pipeline expects the following columns in the csv file: "cluster_id" (contains the cell type names for which differential expression has been computed (same naming scheme as in the config.yaml file)), "gene" (contains the gene names) and "category" (which should be either "Sig. upregulated",  "Sig. downregulated" or "Not deregulated" ).
After editing the config.yaml file the pipeline can be executed with:

~~~
snakemake --use-conda --cores <Ncores> > snakemake.log 2>&1
~~~

