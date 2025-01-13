from re import L
import sys
import os
#os.system("nvidia-smi")
os.environ["NUMBA_CACHE_DIR"] = os.getcwd()
import numba
import numpy
import scvi
import scanpy as sc
import torch
from anndata import read_h5ad,AnnData
from scvi.model.utils import mde
import matplotlib.pyplot as plt
#import pymde
import numpy as np
import pandas as pd

import argparse
# setup some Docker image including argparse
import sys
#import pymde
# maflotho/scvi-tools-gpu:pytorch
#ew
#

def perform_scvi_scanvi(adata, organism,batch_key, labels_key ,unlabeled_category,out_folder, max_epochs,scanvi,n_layers=2, n_latent=30,gene_likelihood="nb"):
    """perform scvi_scanvi

    Args:
        adata (anndata): anndata object
        organism (string): human or mouse
        batch_key (string): batch information obs key
        labels_key (string): cell label information obs key (only for scANVI)
        unlabeled_category (string): Unlabelled category for unknown cell types (only for scANVI)
        out_folder (string): output directory
        max_epochs (int): max number of epochs
        scanvi (bool): perform scANVI after scVI?
        n_layers (int, optional): see scVI documentation. Defaults to 2.
        n_latent (int, optional): see scVI documentation. Defaults to 30.
        gene_likelihood (str, optional): see scVI documentation. Defaults to "nb".
    """
    # setup data for scVI and perform training, parameter tuning might be still necessary
    #batchsize = 1028
    scvi.settings.dl_num_workers = 100
    scvi.settings.num_threads = 100
    
    adata.raw = adata.copy()
    adata.X = adata.X.tocsr()
    adata.layers["counts"] = adata.layers["counts"].tocsr()
    
    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key, layer="counts")
    
    vae = scvi.model.SCVI(adata, n_layers=n_layers, n_latent=n_latent, gene_likelihood=gene_likelihood)
    vae.train(max_epochs=max_epochs,early_stopping = True,use_gpu = True)#batch_size = batchsize
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    adata.layers["scvi_counts"] = vae.get_normalized_expression()
    adata.write_h5ad(f"{out_folder}/scVI_output_{organism}_{batch_key}_{max_epochs}.h5ad")
    vae.save("scVI_model")
    adata.layers["scVI_normalized"] = vae.get_normalized_expression(n_samples=10)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)
    sc.tl.umap(adata)
    adata.obsm['X_umap_raw'] = adata.obsm['X_umap'].copy()
    sc.pp.neighbors(adata, use_rep="X_scVI",key_added = "x_scVI_neighbors")
    sc.tl.leiden(adata,neighbors_key = "x_scVI_neighbors")
    sc.tl.umap(adata,neighbors_key = "x_scVI_neighbors")
    adata.obsm['X_umap_scVI'] = adata.obsm['X_umap'].copy()

    adata.write_h5ad(f"scVI_output_{organism}_{batch_key}_{max_epochs}.h5ad")

    #save the modelled adata
    #start scANVI
    if not scanvi:
        return
    scvi.model.SCANVI.setup_anndata(adata, labels_key=labels_key,unlabeled_category = unlabeled_category,batch_key=batch_key, layer="counts")
    lvae = scvi.model.SCANVI.from_scvi_model(vae,adata=adata,labels_key=labels_key,unlabeled_category=unlabeled_category)   
    lvae.train(max_epochs=max_epochs,early_stopping = True)

    adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
    adata.obs["scanvi_prediction"] = lvae.predict()

    adata.write_h5ad(f"{out_folder}/scANVI_output_{organism}_{batch_key}_{max_epochs}.h5ad")



if __name__ == "__main__":
    

    parser = argparse.ArgumentParser()
    parser.add_argument("file",help = "input file")
    parser.add_argument("id", help = "output id")
    parser.add_argument("--output_folder", help="where to store the files?",default = "./")
    parser.add_argument("--unid_key", help="which class is undefined for scANVI",default = "unID")
    parser.add_argument("--labels_key", help="which class label should be predicted using scANVI",default = "cell_type1")
    parser.add_argument("--batch_key", help="which batch_key should be used by scANVI and scVI",default="batch")
    parser.add_argument("--max_epochs", help="how many epochs should be run?",default=1000)
    parser.add_argument("--scANVI", help="how many epochs should be run?",default=True)

    args = parser.parse_args()

    print(args.file)
    print(args.id)

    if args.unid_key:
        print(args.unid_key)
        
    if not torch.cuda.is_available():
        print("cuda doesn't work, shutting down")
        exit()

    adata = read_h5ad(args.file)
    print(adata)

    perform_scvi_scanvi(adata, organism = args.id, 
                        batch_key= args.batch_key ,
                        max_epochs= args.max_epochs,
                        labels_key= args.labels_key,
                        unlabeled_category = args.unid_key,
                        out_folder = args.output_folder,
                        scanvi = args.scANVI)
