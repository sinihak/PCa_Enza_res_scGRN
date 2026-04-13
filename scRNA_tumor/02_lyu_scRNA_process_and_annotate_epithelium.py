import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
import scipy.io as sio
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import utils as ut
import scipy.cluster.hierarchy as sch
import h5py
import os
import random
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text

random.seed(50932)


from pathlib import Path
datadir = Path.cwd()

adata_concat = ad.read_h5ad(f"{datadir}/scRNA_tumor/anndata/annotated_adata_obj_147180_cells_250813.h5ad")
adata_concat.obs["type"] = adata_concat.obs["type"].cat.reorder_categories(["PC","HSPC","CRPC"])
adata_concat.obs["tumortype_celltype"] = adata_concat.obs["type"].astype("str") + "_" + adata_concat.obs["cell_type"].astype("str")

# Most CPRC epithelial cells in this dataset are from one sample so we are only focusing on PCa and mHSPC

adata_concat = adata_concat[adata_concat.obs["type"].isin(["PC", "HSPC"])]

adata_concat.obs_names = [
    f"{sample}_{cell}" 
    for sample, cell in zip(adata_concat.obs["sample"], adata_concat.obs_names)
]

sio.mmwrite(f"{datadir}/scRNA_tumor/outputs/lyu_scRNA_geneMatrix.mtx",adata_concat.layers["log1p"])
pd.DataFrame(adata_concat.var_names).to_csv(f"{datadir}/scRNA_tumor/outputs/lyu_scRNA_features.csv")
pd.DataFrame(adata_concat.obs_names).to_csv(f"{datadir}/scRNA_tumor/outputs/lyu_scRNA_cells.csv")
pd.DataFrame(adata_concat.obs).to_csv(f"{datadir}/scRNA_tumor/outputs/lyu_scRNA_metadata.csv")


# keep the epithelial cells
adata = adata_concat[adata_concat.obs["cell_type"] == "Epithelial cell"]

# fine-tune the annotation of epithelial cell types:
fig, ax = plt.subplots()
sc.pl.embedding(adata, basis="X_umap", color="leiden", show=False, ax=ax, size=3)
plt.tight_layout()
plt.savefig(f"{datadir}/scRNA_tumor/figures/scRNA_Lyu_umap_epithelium_leiden.pdf")
plt.close()


# cell type-specific marker genes (Supplementary Table 4)
markers = pd.read_csv(f"{datadir}/scRNA_tumor/data/cell_markers_for_plotting.csv", sep=";")

present_genes = set(adata.var_names)  
annot = markers.groupby("cell_type_toPlot")["gene_id_hgnc"].apply(list).to_dict()
epithelial_cells = ["Club", "Basal", "Hillock", "Luminal"]
annot = {key: annot[key] for key in epithelial_cells}


marker_genes_filtered = {
    cell_type: [gene for gene in genes if gene in present_genes]
    for cell_type, genes in annot.items()
}

# Remove cell types that have no remaining genes after filtering
markers= {k: v for k, v in marker_genes_filtered.items() if v}

dp1 = sc.pl.DotPlot(adata,var_names=annot, groupby="leiden",standard_scale="var")
dp1.legend(colorbar_title="Mean expression")
dp1.style(cmap="magma", smallest_dot=30, largest_dot=300)
dp1.savefig(f"{datadir}/scRNA_tumor/figures/scRNA_Lyu_epithelial_markers_dotplot.pdf",bbox_inches="tight")


celltypes = adata.obs["leiden"].astype(str).copy()

celltypes = celltypes.replace({"14": "basal", "5": "club-like"})
celltypes = celltypes.mask(~celltypes.isin(["basal", "club-like"]), "luminal-like")
# convert back to categorical and assign to the anndata object
celltypes = celltypes.astype("category")
adata.obs["cell_type_refined"] = celltypes

adata.obs["tumortype_celltype"] = adata.obs["type"].astype("str") + "_" + adata.obs["cell_type_refined"].astype("str")

# map the refined epithelial cell types back to the full data
adata_concat.obs["cell_type_refined"] = adata_concat.obs["cell_type"].astype(str)
adata_concat.obs.loc[adata.obs_names, "cell_type_refined"] = adata.obs["cell_type_refined"]
adata_concat.obs["cell_type_refined"] = adata_concat.obs["cell_type_refined"].astype("category")

adata_concat.obs["tumortype_celltype"] = adata_concat.obs["type"].astype("str") + "_" + adata_concat.obs["cell_type_refined"].astype("str")
adata_concat.obs["tumortype_celltype"] = adata_concat.obs["tumortype_celltype"].astype("category")

adata_concat.write(f"{datadir}/scRNA_tumor/anndata/annotated_adata_obj_refined.h5ad", compression="gzip")


cats = adata_concat.obs["tumortype_celltype"].cat.categories.tolist()
# Group by prefix (case-insensitive, safer)
pca  = [c for c in cats if c.startswith("PC")]
hspc = [c for c in cats if c.startswith("HSPC")]

# Combine in desired order
order = pca + hspc


# ----------------------------------------------
# Visualize SOX4 regulon (supplementary Fig 2I) 
# ----------------------------------------------

sox4_regulon = ["SESN3","CD24","KCNMB4","MBNL2","ZMAT1","TUBA1A","MARCKS","ZSWIM5","TMSB4X","C16orf89","ZNF608","MET","SOX4"]
X_sox4 = adata[:, sox4_regulon].X.toarray()
sox4_linkage = sch.linkage(X_sox4.T, method="ward")  # Cluster columns
sox4_dendro = sch.dendrogram(sox4_linkage, no_plot=True)["leaves"]
sox4_order = [sox4_regulon[i] for i in sox4_dendro]

dp2 = sc.pl.DotPlot(adata,var_names=sox4_order, groupby="tumortype_celltype",standard_scale="var",categories_order=order)
dp2.legend(colorbar_title="Mean expression")
dp2.style(cmap="magma", smallest_dot=30, largest_dot=300)
dp2.savefig(f"{datadir}/scRNA_tumor/figures/scRNA_Lyu_tumor_type_refined_epithelial_cell_type_SOX4_regulon_dotplot.pdf",bbox_inches="tight", dpi=300)
