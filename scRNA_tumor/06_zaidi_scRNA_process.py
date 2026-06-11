
import anndata as ad
import pandas as pd
import scanpy as sc
import seaborn as sns
import doubletdetection
import numpy as np
import scipy.io as sio
import random
from scipy.io import mmread
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from pathlib import Path
datadir = Path.cwd() 

random.seed(129412)

# In R 
# library(Matrix)
# seu <- readRDS(paste0(datadir, "/scRNA_tumor/rds/GSE264573_msk.integrated.remove.cellcycle.allcells.rds"))


# counts <- seu[['RNA']]$counts
# cells <- colnames(counts)
# genes <- rownames(counts)

# writeMM(counts, paste0(datadir,"/scRNA_tumor/outputs/rna_counts.mtx"))
# write.csv(cells, paste0(datadir,"/scRNA_tumor/outputs/rna_cells.csv"), quote =FALSE)
# write.csv(genes, paste0(datadir,"/scRNA_tumor/outputs/rna_genes.csv"), quote= FALSE)
# write.csv(seu@meta.data, paste0(datadir,"/scRNA_tumor/outputs/rna_metadata.csv"), quote=FALSE)

# counts = mmread(f"{datadir}/scRNA_tumor/outputs/rna_counts.mtx")


# Cell and feature information
# cells = pd.read_csv(f"{datadir}/scRNA_tumor/outputs/rna_cells.csv", index_col=0).iloc[:, 0]
# features = pd.read_csv(f"{datadir}/scRNA_tumor/outputs/rna_genes.csv", index_col=0).iloc[:, 0]

# # Make AnnData object 
# adata = sc.AnnData(counts.T)
# adata.obs_names = cells
# adata.var_names = features


# # Add metadata
# cell_meta = pd.read_csv(f"{datadir}/scRNA_tumor/outputs/rna_metadata.csv", index_col=0).loc[adata.obs_names, : ]

# for col in cell_meta.columns:
#     adata.obs[col] = cell_meta[col].values

# adata.X = adata.X.tocsr()

# adata.layers["counts"] = adata.X.copy()

# adata.write(f"{datadir}/anndata/zaidi_scRNA_anndata.h5ad", compression='gzip')

adata = ad.read_h5ad(f"{datadir}/scRNA_tumor/anndata/zaidi_scRNA_anndata.h5ad")


plt.figure(figsize=(4,4))
sc.pl.violin(adata,["nCount_RNA","nFeature_RNA" ,"percent.mt"],
    jitter=0.4,
    multi_panel=True,
    stripplot=False,
)
plt.savefig(f"{datadir}/scRNA_tumor/figures/zaidi_rna_QC_violins_prefilt.pdf", dpi=300, bbox_inches="tight", pad_inches=0.5)

sc.pp.filter_genes(adata, min_cells=500)

adata = adata[(adata.obs.nCount_RNA < 80000) & (adata.obs.nCount_RNA >= 1000)  & (adata.obs.nFeature_RNA < 12000),:,].copy()

# adata.var_names[adata.var_names.str.startswith("MT-")]

plt.figure(figsize=(4,4))
sc.pl.violin(adata,["nCount_RNA","nFeature_RNA"],
    jitter=0.4,
    multi_panel=True,
    stripplot=False,
)
plt.savefig(f"{datadir}/scRNA_tumor/figures/zaidi_rna_QC_violins.pdf", dpi=300, bbox_inches="tight", pad_inches=0.5)

# perform doublet detection
clf = doubletdetection.BoostClassifier(
    n_iters=10,
    clustering_algorithm="louvain",
    standard_scaling=True,
    pseudocount=0.1,
    n_jobs=-1,
)

doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
doublet_score = clf.doublet_score()

p1 = doubletdetection.plot.convergence(clf, show=True, p_thresh=1e-16, voter_thresh=0.5)
p1.savefig(f'{datadir}/scRNA_tumor/figures/zaidi_doublet_detection_convergence_test.pdf')

adata.obs["doublet"] = doublets
adata.obs["doublet_score"] = doublet_score

print(adata.obs['doublet'].value_counts())

adata = adata[adata.obs['doublet'] != 1.0]

adata.write(f"{datadir}/scRNA_tumor/anndata/zaidi_scRNA_anndata_filtered.h5ad", compression='gzip')

sc.tl.pca(adata, n_comps=100, svd_solver="auto")

# subset to epithelial only
adata = adata[adata.obs.coarse_ano == 'Epi_Neuroendo'].copy()

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")

sc.pl.highly_variable_genes(adata, log=True)
plt.savefig(f"{datadir}/scRNA_tumor/figures/zaidi_hvg_plot.pdf")
plt.tight_layout()
plt.close()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.write(f"{datadir}/scRNA_tumor/anndata/zaidi_scRNA_anndata_unscaled.h5ad", compression='gzip')
adata = ad.read_h5ad(f"{datadir}/scRNA_tumor/anndata/zaidi_scRNA_anndata_unscaled.h5ad")

sc.pp.scale(adata)
sc.tl.pca(adata, n_comps=100, svd_solver="auto")

sc.pl.pca_variance_ratio(adata, show=False, n_pcs=50)
plt.savefig(f"{datadir}/figures/zaidi_pca_rankplot.pdf")
plt.close()


adata.write(f"{datadir}/scRNA_tumor/anndata/zaidi_scRNA_anndata_normalized.h5ad", compression='gzip')
adata = ad.read_h5ad(f"{datadir}/scRNA_tumor/anndata/zaidi_scRNA_anndata_normalized.h5ad")

sc.pp.neighbors(adata, metric="cosine",n_pcs=20, n_neighbors=50)

sc.tl.umap(adata)

to_plot = ['sample', 'subtype', 'coarse_ano', 'patient_prev', 'nCount_RNA', 'percent.mt','nFeature_RNA']

filename = (f"{datadir}/scRNA_tumor/figures/zaidi_rna_umaps.pdf")

with PdfPages(filename) as pdf:
    for feature in to_plot:
        sc.pl.umap(adata, size=10, color=feature, show=False)
        plt.tight_layout()
        plt.title(feature)
        pdf.savefig(dpi=300,bbox_inches="tight", pad_inches=0.5)
        plt.close()


resolutions = [0.1,0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
n_clusters = []
for res in resolutions:
    sc.tl.leiden(adata, resolution=res)
    n_clusters.append(adata.obs['leiden'].nunique())

plt.figure(figsize=(4,4))
plt.plot(resolutions, n_clusters, marker='o')
plt.xlabel('Resolution')
plt.ylabel('Number of clusters')
plt.tight_layout()
plt.savefig(f'{datadir}/scRNA_tumor/figures/zaidi_rna_clustering_resolution.pdf')
plt.close()

sc.tl.leiden(adata, resolution=0.6)

sc.pl.pca(adata, color='leiden', size=10)
plt.tight_layout()
plt.gcf().savefig(f"{datadir}/scRNA_tumor/figures/zaidi_rna_pca_leiden.pdf")
plt.close()

adata.write(f"{datadir}/scRNA_tumor/anndata/zaidi_scRNA_anndata_clustered.h5ad", compression='gzip')

# for AUCEll
# sio.mmwrite(f"{datadir}/outputs/zaidi_scRNA_geneMatrix_epithelium.mtx",adata.layers["counts"])
# pd.DataFrame(adata.var_names).to_csv(f"{datadir}/outputs/zaidi_scRNA_features_epithelium.csv")
# pd.DataFrame(adata.obs_names).to_csv(f"{datadir}/outputs/zaidi_scRNA_cells_epithelium.csv")
# pd.DataFrame(adata.obs).to_csv(f"{datadir}/outputs/zaidi_scRNA_metadata_epithelium.csv")
