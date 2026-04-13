import scglue
import re
import anndata
import torch
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import rcParams
from pathlib import Path
datadir = Path.cwd() 


SEED = 1234
lr=0.0002
date = "20240307" 

rna = anndata.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/rna_{date}.h5ad")
atac = anndata.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/atac_{date}.h5ad")

# Read the guidance graph
guidance = nx.read_graphml(f"{datadir}/TF_regulatory_inference/outputs/dist_{date}.graphml.gz")

print("Configuring data...")
scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca",
)

hvg_reachable = scglue.graph.reachable_vertices(guidance, rna.var.query("highly_variable").index)
atac.var["highly_variable"] = [item in hvg_reachable for item in atac.var_names]

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
    )

atac.var["highly_variable"].sum()

print("Subsetting graph for HVGs...")
from itertools import chain
guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()

print("Training the model...")
glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, guidance_hvf,
    fit_kws={"directory":outdir + "glue-"+date,
             "safe_burnin":False,"align_burnin":np.inf},
    init_kws={"latent_dim":50},
    compile_kws={"lr":lr}
)

print("Saving the model...")
glue.save(f"{datadir}/TF_regulatory_inference/outputs/glue_{date}.dill")
