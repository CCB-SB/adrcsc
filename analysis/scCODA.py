import importlib
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
import tensorflow as tf
import sccoda.datasets as scd

ct_ref = "CD4+ T cell"

cell_counts = pd.read_csv('scCODA/CellCounts_Apoe_broad.csv')

data_salm = dat.from_pandas(cell_counts, covariate_columns=["SCMD", "Age", "ApoE", "Sex", "Diagnosis"])
model_salm = mod.CompositionalAnalysis(data_salm, formula="C(Diagnosis, Treatment('HC')) + Age + ApoE + Sex", reference_cell_type=ct_ref)
sim_results = model_salm.sample_hmc()
sim_results.set_fdr(est_fdr=0.05)
res  = sim_results.summary()
effects_res = sim_results.effect_df
path = "scCODA/broadAnnotation/Output_All.csv"
effects_res.to_csv(path)
f = open("scCODA/broadAnnotation/All.txt", "a")
f.write(pd.Series.to_string(sim_results.credible_effects()))
f.close()


data_salm = dat.from_pandas(cell_counts, covariate_columns=["SCMD", "Age", "ApoE", "Sex", "Diagnosis"])
data_salm = data_salm[data_salm.obs["Sex"].isin(["male"])]
model_salm = mod.CompositionalAnalysis(data_salm, formula="C(Diagnosis, Treatment('HC')) + Age + ApoE", reference_cell_type=ct_ref)
sim_results = model_salm.sample_hmc()
sim_results.set_fdr(est_fdr=0.05)
res  = sim_results.summary()
effects_res = sim_results.effect_df
path = "scCODA/broadAnnotation/Output_Male.csv"
effects_res.to_csv(path)
f = open("scCODA/broadAnnotation/Male.txt", "a")
f.write(pd.Series.to_string(sim_results.credible_effects()))
f.close()


data_salm = dat.from_pandas(cell_counts, covariate_columns=["SCMD", "Age", "ApoE", "Sex", "Diagnosis"])
data_salm = data_salm[data_salm.obs["Sex"].isin(["female"])]
model_salm = mod.CompositionalAnalysis(data_salm, formula="C(Diagnosis, Treatment('HC')) + Age + ApoE", reference_cell_type=ct_ref)
sim_results = model_salm.sample_hmc()
sim_results.set_fdr(est_fdr=0.05)
res  = sim_results.summary()
effects_res = sim_results.effect_df
path = "scCODA/broadAnnotation/Output_Female.csv"
effects_res.to_csv(path)
f = open("scCODA/broadAnnotation/Female.txt", "a")
f.write(pd.Series.to_string(sim_results.credible_effects()))
f.close()



ct_ref = "Treg CD4+ cell"
cell_counts = pd.read_csv('scCODA/CellCounts_Apoe_fine.csv')

data_salm = dat.from_pandas(cell_counts, covariate_columns=["SCMD", "Age", "ApoE", "Sex", "Diagnosis"])
model_salm = mod.CompositionalAnalysis(data_salm, formula="C(Diagnosis, Treatment('HC')) + Age + ApoE + Sex", reference_cell_type=ct_ref)
sim_results = model_salm.sample_hmc()
sim_results.set_fdr(est_fdr=0.05)
res  = sim_results.summary()
effects_res = sim_results.effect_df
path = "scCODA/fineAnnotation/Output_All.csv"
effects_res.to_csv(path)
f = open("scCODA/fineAnnotation/All.txt", "a")
f.write(pd.Series.to_string(sim_results.credible_effects()))
f.close()


data_salm = dat.from_pandas(cell_counts, covariate_columns=["SCMD", "Age", "ApoE", "Sex", "Diagnosis"])
data_salm = data_salm[data_salm.obs["Sex"].isin(["male"])]
model_salm = mod.CompositionalAnalysis(data_salm, formula="C(Diagnosis, Treatment('HC')) + Age + ApoE", reference_cell_type=ct_ref)
sim_results = model_salm.sample_hmc()
sim_results.set_fdr(est_fdr=0.05)
res  = sim_results.summary()
effects_res = sim_results.effect_df
path = "scCODA/fineAnnotation/Output_Male.csv"
effects_res.to_csv(path)
f = open("scCODA/fineAnnotation/Male.txt", "a")
f.write(pd.Series.to_string(sim_results.credible_effects()))
f.close()


data_salm = dat.from_pandas(cell_counts, covariate_columns=["SCMD", "Age", "ApoE", "Sex", "Diagnosis"])
data_salm = data_salm[data_salm.obs["Sex"].isin(["female"])]
model_salm = mod.CompositionalAnalysis(data_salm, formula="C(Diagnosis, Treatment('HC')) + Age + ApoE", reference_cell_type=ct_ref)
sim_results = model_salm.sample_hmc()
sim_results.set_fdr(est_fdr=0.05)
res  = sim_results.summary()
effects_res = sim_results.effect_df
path = "scCODA/fineAnnotation/Output_Female.csv"
effects_res.to_csv(path)
f = open("scCODA/fineAnnotation/Female.txt", "a")
f.write(pd.Series.to_string(sim_results.credible_effects()))
f.close()
