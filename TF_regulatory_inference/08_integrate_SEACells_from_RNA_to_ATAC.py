import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd
import scglue 
import seaborn as sns
import scipy
from IPython import display
from matplotlib import rcParams
from networkx.algorithms.bipartite import biadjacency_matrix
from networkx.drawing.nx_agraph import graphviz_layout
import scanpy as sc
import matplotlib.pyplot as plt
import pylab as pl
import random

from pathlib import Path
# fetches current wd
datadir = Path.cwd() 

random.seed(321)
date= "20240307"
num="250" # nCells/metacell

scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)

rna = ad.read_h5ad(os.path.join(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-w-leiden-phases-{date}.h5ad"))
atac = ad.read_h5ad(os.path.join(f"{datadir}/TF_regulatory_inference/anndata/atac-emb-w-leiden-{date}.h5ad"))

combined = ad.read_h5ad(os.path.join(f"{datadir}/TF_regulatory_inference/anndata/combined-emb-{date}.h5ad"))

sc_rna = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_rna_seacells_mod_num_{num}_{date}.csv"))

sc_rna.set_index("index", inplace=True)
sc_rna = sc_rna.reindex(rna.obs_names)


rna.obs["SEACell"] = sc_rna["SEACell"]

# transfer "SEAcell" from atac to rna using GLUE embedding
scglue.data.transfer_labels(rna, atac, "SEACell", use_rep="X_glue")

df = pd.DataFrame(data=atac.obs)

modalities = ["RNA", "ATAC"]
# Create a dictionary to store the split AnnData objects
ads = {}
# Iterate over each specified modality and create a subset AnnData
for modality in modalities:
    combined_sub = combined[combined.obs["modality"] == modality].copy()
    ads[modality] = combined_sub

rna_comb = ads["RNA"]
atac_comb = ads["ATAC"]

rna_comb.obs["SEACell"] = rna.obs["SEACell"]
atac_comb.obs["SEACell"] = atac.obs["SEACell"]


combined_back = ad.concat([rna_comb, atac_comb], axis=0, join="outer", label="modality", keys=["RNA", "ATAC"])

sc_atac = atac.obs["SEACell"]
print(len(sc_atac.unique()))
print(len(sc_rna.unique()))

# write seacells into csv file

sc_atac.to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_atac_seacells_integrated_num_{num}_{date}.csv")
atac.write(os.path.join(seadir,f"/TF_regulatory_inference/anndata/lncap_atac_seacells_integrated_num_{num}_{date}.h5ad"))

sc_rna.to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_rna_seacells_integrated_num_{num}_{date}.csv")
rna.write(f"{datadir}/TF_regulatory_inference/anndata/lncap_rna_seacells_integrated_num_{num}_{date}.h5ad")

rna.obs.to_csv(f"{datadir}/TF_regulatory_inference/lncap_rna_seacells_metadata_integrated_num_{num}_{date}.csv")




