import gc
import time
from sc_atlas import *
import scanpy as sc
from scanpy import pp, tl, pl
import json
import argparse


"""Preprocessing 

    1. doublet removal
    2. merging
    3. quality control


    input : 
        - raw_atlas.h5
        - organism
        - protocol
        - excluded datasets
"""
def preprocessingI_from_scratch(atlas, processed_root, out, method, organism, from_scratch = False,exclude = []):
    print("Excluding :",exclude)
    atlas_processed = []
    if from_scratch:
        atlas_processed = ScAtlas(file = processed_root, method_lvl = ["10x","10x_sn"],organism_lvl=["human","mouse"],mode = "w",overwrite = True)
    else:
        atlas_processed = ScAtlas(file = processed_root, mode = "r+")
        
    atlas.perform_preprocessingI(out_atlas = atlas_processed,method = method, organism = organism,exclude = exclude)
    atlas.export_scrubbed_matrix(atlas_processed,method = method, organism = organism,exclude = exclude)
    atlas.f.close()
    gc.collect()
    atlas_processed.merge(method,organism,out =out,exclude = exclude)
    atlas_processed.get_h5ad(out, file = f"{out}_{organism[0]}.h5ad")


def preprocessingII(file,out,mt = 'mt',dim_reduction = True, exclude =[]):
    """ remove less abundant genes, mitochondrial genes etc

    Args:
        file (string): _description_
        out (string): _description_
        mt (str, optional): _description_. Defaults to 'mt'.
    """   
    print(file)
    adata = sc.read_h5ad(file)
    for ds in exclude:
        adata = adata[adata.obs["dataset"] != ds]
    file  = file.replace("/","_")
    #adata.var_names_make_unique() 
    #adata.obs_names_make_unique() 
    cortex,noncortex = adata[adata.obs["tissue"] == "Cortex"], adata[adata.obs["tissue"] != "Cortex"]
    print(np.unique(adata.obs["tissue"]))
    print(np.unique(adata.obs["cell_type1"]))

    print(cortex)

    def process(adata, out):
        sc.pl.highest_expr_genes(adata, n_top=20, save= f"{file}.png")
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        adata.var["mt"] = adata.var_names.str.startswith(f'{mt}-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                  multi_panel=True,jitter = False, save = f"{file}_violin.png",)
        sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',save = f"{file}_pct_counts.png")
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',save= f"{file}_genes_by_counts.png")
        #adata.raw = adata
        adata.layers["counts"] = adata.X.copy()
        adata = adata[adata.obs.n_genes_by_counts < 7500, :]
        adata = adata[adata.obs.pct_counts_mt < 5, :]
        adata.write_h5ad(f"{out}_raw.h5ad")
        #sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pl.highly_variable_genes(adata,save = f"{file}_high_var.png")
        adata = adata[:, adata.var.highly_variable]
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        AnnData.write_h5ad(adata,f"{out}.h5ad")

    process(noncortex,f"{out}_noncortex")
    process(cortex,f"{out}_cortex")



def run(prep,organism,method,root,from_scratch = False):
    
    """ run Perform pre-processing 1 and 2:
            1. mapping and appending matrices + perform scrublet per study
            2. filter low quality cells and split by cortex/non-cortex region
    Args:
        prep (list(bool)): which processing to perform
        organism (string): "human" or "mouse"
        method (string): sequencing protocol "10x" or "10x_sn"
        root (string) : path to atlas h5 file
        from_scratch (bool) : replace intermediate files
        ": _description_. Defaults to 'mt'.
    """  
    
    p_root = f"{root}atlas_processed.h5"
    id = f"{organism[0]}_{method[0]}"

    outI = f"merged_{id}"
    outII = f"merged_{id}.h5ad"

    prepI, prepII = prep
    exclude = {"grubman_2019","zheng_2017","tabula_muris_senis"}
    
    mt = "mt"
    if organism == ["human"]: mt = "MT"

    if prepI :
        atlas = ScAtlas(root,mode ="r+")
        print("perform preprocessing I")
        preprocessingI_from_scratch(atlas,processed_root = p_root,out=outI,method=method,organism=organism,from_scratch = from_scratch,exclude = exclude)
        gc.collect()
    if prepII:
        print("perform preprocessing II")
        preprocessingII(f"{outI}.h5ad",out = outII,mt = mt,exclude = exclude)
        gc.collect()

    print(f"Took :\t{(time.time() - start)/60.0} minutes")



if __name__ == "__main__":
    start = time.time()
    
    parser = argparse.ArgumentParser()    
    parser.add_argument("atlas", help="h5 atlas file")
    parser.add_argument("--prepI", help="perform preprocessing I?",action='store_true')
    parser.add_argument("--prepII", help="perform preprocessing II?",action='store_true')
    parser.add_argument("--method",help = "10x or 10x_sn",default = "10x")
    parser.add_argument("--organism",help = "10x or 10x_sn",default = "mouse")
    parser.add_argument("--from_scratch",help = "10x or 10x_sn",action='store_true')
    parser.add_argument("--config", help="dpi",default = "./config_files/core3_config.json")
    args = parser.parse_args()
    config = json.load(open(args.config))

    args = parser.parse_args()
    method =[args.method]
    organism = [args.organism]
    run([args.prepI,args.prepII],organism=organism,method=method, root = config["root"], from_scratch = args.from_scratch)    


