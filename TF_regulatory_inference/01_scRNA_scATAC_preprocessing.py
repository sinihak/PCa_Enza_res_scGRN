import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
import pandas as pd
import numpy as np
from matplotlib import rcParams
from typing import Optional
from itertools import chain
from scipy.io import mmread
import itertools
import seaborn as sns
from tqdm import tqdm
import networkx as nx
from networkx.algorithms.bipartite import biadjacency_matrix
from pathlib import Path

datadir = Path.cwd() 

scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)

date="20240307" 
print("processing RNA data")
hvgs=4000

# read the cells to be removed:
toRM = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/scRNA_cell_barcodes_to_remove.csv")
# Extract the cell names from the DataFrame
toRM = toRM["cell_names"].tolist()

countsRNA = mmread(f"{datadir}/TF_regulatory_inference/outputs/rna_counts_{date}.mtx")
# Cell and feature information
cellsRNA = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/cell_rna_{date}.csv", index_col=0).iloc[:, 0]
features = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/genes_rna_{date}.csv", index_col=0).iloc[:, 0]

# Make AnnData object (RNA)
rna = sc.AnnData(countsRNA.T)
rna.obs_names = cellsRNA
rna.var_names = features

# Add metadata
cell_meta = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/metadata_rna_{date}.csv", index_col=0).loc[rna.obs_names, : ]
for col in cell_meta.columns:
    rna.obs[col] = cell_meta[col].values

rna.X = rna.X.tocsr()

rna.layers["counts"] = rna.X.copy()

print(len(rna.obs_names))
# remove the specific cells:
rna = rna[~rna.obs.index.isin(toRM)].copy()
print(len(rna.obs_names))

## Add promoter annotation to RNA
scglue.data.get_gene_annotation(
    rna, gtf=f"{datadir}/data/gencode.v41.chr_patch_hapl_scaff.annotation.gtf.gz",
    gtf_by="gene_name"
)

print(rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head(10))

# Remove NA cases and make start/end to integer
keepGenes = rna.var.index[rna.var["chrom"].notna()]
keep = np.in1d(rna.var_names, keepGenes)
rna = rna[:,keep]
rna.var.loc[:,["chromStart","chromEnd"]] = rna.var.loc[:,["chromStart","chromEnd"]].fillna(0).astype(int)
print(rna.var.head(10))

# Get HVG
sc.pp.highly_variable_genes(rna, n_top_genes=hvgs, flavor="seurat_v3")

print(rna.var["highly_variable"].sum())

# Normalisation, scaling
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)

# Dim reduction
sc.tl.pca(rna, n_comps=50, svd_solver="auto")
sc.pp.neighbors(rna, metric="cosine")
sc.tl.umap(rna)
# rna.obs["seurat_clusters"] = rna.obs["seurat_clusters"].astype(str)
# sc.pl.umap(rna, color=["seurat_clusters"], save="_lncap_scRNA_UMAP_{date}.pdf")

print("processing ATAC data")

# read cell and peak information

countsATAC = mmread(f"{datadir}/TF_regulatory_inference/outputs/atac_counts_{date}.mtx")
cellsATAC = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/cell_atac_{date}.csv"), index_col=0.iloc[:, 0]
peaks = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/peaks_atac_{date}.csv"), index_col=0.iloc[:, 0]

# Make AnnData object (ATAC)
atac = sc.AnnData(countsATAC.T)
atac.obs_names = cellsATAC
atac.var_names = peaks

# Add metadata
cell_meta = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/metadata_atac_{date}.csv"), index_col=0).loc[atac.obs_names, : ]
for col in cell_meta.columns:
    atac.obs[col] = cell_meta[col].values

atac.X = atac.X.tocsr()

# lsi/neighbors/umap
scglue.data.lsi(atac, n_components=50, n_iter=15)
sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
sc.tl.umap(atac)
# atac.obs["seurat_clusters"] = atac.obs["seurat_clusters"].astype(str)
# sc.pl.umap(atac,color=["seurat_clusters"], save="_lncap_scATAC_UMAP_{date}.pdf")

# format atac
split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
atac.var.head()

# get genes, peaks, tss, and promoters
genes = scglue.genomics.Bed(rna.var.assign(name=rna.var_names))
peaks = scglue.genomics.Bed(atac.var.assign(name=atac.var_names))
tss = genes.strand_specific_start_site()
promoters = tss.expand(2000, 0)

# Build genomic distance graph with power decay law weighting (+/- 150kb)
dist_graph = scglue.genomics.window_graph(
    promoters, peaks, 150000,
    attr_fn=lambda l, r, d: {
        "dist": abs(d),
        "weight": scglue.genomics.dist_power_decay(abs(d)),
        "type": "dist"
    }
)

dist_graph = nx.DiGraph(dist_graph)
dist_graph.number_of_edges()

hvg_reachable = scglue.graph.reachable_vertices(dist_graph, rna.var.query("highly_variable").index)

# modify highly variable in ATAC
atac.var["highly_variable"] = [item in hvg_reachable for item in atac.var_names]
print(atac.var["highly_variable"].sum())

dist_graph = scglue.graph.compose_multigraph(dist_graph, dist_graph.reverse())

for item in itertools.chain(atac.var_names, rna.var_names):
    dist_graph.add_edge(item, item, weight=1.0, type="self-loop")
nx.set_edge_attributes(dist_graph, 1, "sign")

scglue.graph.check_graph(dist_graph, [rna, atac])

rna.var["artif_dupl"] = rna.var["artif_dupl"].astype(str) # convert to string to avoid error

print(len(rna.obs_names))
print(len(atac.obs_names))
print("saving the output")

rna.write(f"{datadir}/TF_regulatory_inference/anndata/rna_{date}.h5ad")
atac.write(f"{datadir}/TF_regulatory_inference/anndata/atac_{date}.h5ad")
nx.write_graphml(dist_graph, f"{datadir}/TF_regulatory_inference/outputs/dist_{date}.graphml.gz")

