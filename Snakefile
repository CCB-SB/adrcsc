from os.path import join, isfile, dirname, basename
from glob import glob
from collections import defaultdict

#from script.utils import read_categories

import pandas as pd

conf_file = "config.yaml"

if isfile(conf_file):

    configfile: conf_file


metadata = pd.read_csv("data/metadata.csv", sep="\t", index_col="Date Processed at CG")


celltype_annotations = {
    "celltype_cluster": {
        "cDC1","cDC2","CD16 Monocytes",  
        "Megacaryocytes","DN1 B cell","CD14 Monocytes",
        "Th1/Th2/Th17 cell","CD56-Dim, CD16 NK cell","M-Memory B cell",
        "Proliferating Monocytes","CD8+ Tem T cell","DN B cell",
        "CD4+ Proliferating T cell","CD56-Bright NK cell","Naive CD4+ T cell",
        "HSPC / Progenitor TO_DELETE","CD8+ T cell","Treg CD4+ cell",
        "Transitional B cell","Naive B cell","CD4+ Proliferating",
        "C-Memory B cell","Naive CD8+ T cell","CD4+ Tcm/Tscm T cell",
        "CD8+ Tcm T cell","CD4+ Tem cell","CD8+ Tte T cell",
        "Proliferating T cell","CD4+ Tcm T cell","CD4+ Tfh cell",
        "NKT-like cell","pDC cell","Red Blood Cells"
    },
    "celltype": {
        "CD4+ T-Helper Cell",
	"Naive CD4+ T cell",
	"CD4+ Memory T cell",
	"Gamma delta T cell",
	"Mucosal associated invariant T cell",
	"CD8+ Memory T cell",
	"Naive CD8+ T cell",
	"Proliferating CD8+ T cell",
	"Memory B cell",
	"Double negative B cell",
	"Naive B cell",
	"B cell",
	"Transitional B cell",
	"Conventional Dendritic cell",
	"CD16+ Monocyte",
	"Megacaryocytes",
	"Proliferating Monocyte",
	"CD14+ Monocyte",
	"Plasmacytoid Dendritic cell",
	"NKT-like cell",
	"Proliferating CD4+ T cell",
	"Plasmablasts",
	"CD56-Dim NK cell",
	"CD56-Bright NK cell"
    }
}


celltype_annotations["celltype_w_batch_and_sex"] = celltype_annotations["celltype"]

CATEGORIES = config["enrichment"]["categories"]

COMPARISONS = ["MCIvsHC", "ADvsHC", "PDvsHC", "PDMCIvsHC"]

#SETTING_CELL_POP_UP_ANNOT = defaultdict(list)
#SETTING_CELL_POP_DOWN_ANNOT = defaultdict(list)
#SETTING_CELL_POP_ALL_ANNOT = defaultdict(list)
#SETTING_CELL_POP_GSEA_ANNOT = defaultdict(list)
#SETTING_CELL_POP_GSEA_SHAP_ANNOT = defaultdict(list)
#for c in COMPARISONS:
#    for cannot, celltypes in celltype_annotations.items():
#        for ctype in celltypes:
#            SETTING_CELL_POP_UP_ANNOT[cannot].append("{}|{}|{}_up".format(cannot, c, ctype))
#            SETTING_CELL_POP_DOWN_ANNOT[cannot].append("{}|{}|{}_down".format(cannot, c, ctype))
#            SETTING_CELL_POP_ALL_ANNOT[cannot].append("{}|{}|{}_all".format(cannot, c, ctype))
#            SETTING_CELL_POP_GSEA_ANNOT[cannot].append("{}|{}|{}".format(cannot, c, ctype))
#            if c in SHAP_COMPARISONS:
#                SETTING_CELL_POP_GSEA_SHAP_ANNOT[cannot].append("{}|{}|{}".format(cannot, c, ctype))

#SETTING_CELL_POP = ["{}|{}|{}_{}".format(cannot, c, ctype, regdir) for c in COMPARISONS for cannot, celltypes in celltype_annotations.items() for ctype in celltypes for regdir in ["up", "down", "all"]]

#TARGET_DBS = [cat[0] for cat in read_categories(config['genetrail3']['dge_categories'])]
#TARGET_DBS = ["GO_-_Biological_Process", "GO_-_Cellular_Component", "GO_-_Molecular_Function", "KEGG_-_Pathways", "Reactome_-_Pathways", "Pfam_-_Protein_families", "MSigDB_v7_C6_Oncogenic_Signatures", #"MSigDB_v7_C7_Immunologic_Signatures", "MSigDB_v7_H_Hallmark"]

TOOLS = ("limma_voom") #("edgeR", "limma_trend", "limma_voom", "deseq2")

### No pre-Processing ! ####
rule all:
    input:
        expand(join(config["results_dir"], "de", "heatmap", "ds_pb_{tool}_all_comparisons_{celltype}.all.pdf"), celltype=("celltype_cluster"), tool=TOOLS), #, "celltype_w_batch_and_sex", "celltype_cluster", "celltype_cluster_w_batch_and_sex"



include: "rules/diff_exp.smk"






