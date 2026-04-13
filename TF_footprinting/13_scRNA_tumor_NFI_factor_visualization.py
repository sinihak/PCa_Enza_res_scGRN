
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import utils as ut
import scipy.cluster.hierarchy as sch
import os
import random
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
datadir = Path.cwd() 


# -------------------------------------------------------------------
# NFI factor expression in Lyu et al 2024 data (Supplementary Fig 5C)
# -------------------------------------------------------------------


adata_lyu = ad.read_h5ad(f"{datadir}/scRNA_tissue/anndata/lyu_annotated_adata_obj_refined.h5ad")
adata_lyu.obs["tumortype_celltype"] = adata_lyu.obs["type"].astype("str") + "_" + adata_lyu.obs["cell_type"].astype("str")
adata_lyu.obs["tumortype_celltype"] = adata_lyu.obs["tumortype_celltype"].astype("category")
# Keep the PCa and HSPC cells
# adata_lyu = adata_lyu[adata_lyu.obs["type"].isin(["PC", "HSPC"])]

nfi_genes = ["NFIA", "NFIB", "NFIC", "NFIX"]

cats = adata_lyu.obs["tumortype_celltype"].cat.categories.tolist()
# Group by prefix
pca  = [c for c in cats if c.startswith("PC")]
hspc = [c for c in cats if c.startswith("HSPC")]

# Combine in desired order, PCa cell types before HSPC
order = pca + hspc


genes = ["NFIA", "NFIB", "NFIC", "NFIX"]

dp1 = sc.pl.DotPlot(adata_lyu, var_names=genes, groupby="tumortype_celltype", standard_scale="var")
dp1.style(cmap="magma",smallest_dot=50, largest_dot=500)
dp1.swap_axes()
dp1.savefig(f"{datadir}/scRNA_tissue/figures/scRNA_Lyu_dotplot_NFI.pdf",bbox_inches="tight", dpi=300, pad_inches=0.5)
plt.close()

# -------------------------------------------------------------------
# NFI factor expression in Dong et al 2020 data (Supplementary Fig 5D)
# -------------------------------------------------------------------

adata_dong = ad.read_h5ad(f"{datadir}/scRNA_tissue/anndata/scRNA_dong_adata_epithelium_celltypes_annotated_w_neuroendocrine.h5ad")

dp2 = sc.pl.DotPlot(adata_dong,var_names=genes, groupby="celltypes",standard_scale="var")
dp2.legend(colorbar_title="Mean expression")
dp2.style(cmap="magma", smallest_dot=30, largest_dot=300)
dp2.swap_axes()
dp2.savefig(f"{datadir}/scRNA_tissue/figures/scRNA_Dong_dotplot_NFI.pdf", dpi=300, pad_inches=0.5)
plt.close()



