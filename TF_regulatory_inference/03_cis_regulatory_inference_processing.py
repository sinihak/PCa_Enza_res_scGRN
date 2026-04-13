import os
import re
import seaborn as sns
import anndata as ad
import torch
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import rcParams
import seaborn as sns
import scglue
import matplotlib.pyplot as plt

from pathlib import Path
datadir = Path.cwd() 


SEED = 1234
date = "20240307"
th = 0.1

# Read data
rna = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/rna_{date}.h5ad")
atac = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/atac_{date}.h5ad")
guidance = nx.read_graphml(f"{datadir}/TF_regulatory_inference/outputs/dist_{date}.graphml.gz")

# read the trained model
glue = scglue.models.load_model(f"{datadir}/TF_regulatory_inference/outputs/glue_{date}.dill")

# configure the model
scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)

# fetch the peaks reachable from HVGs and assign them as highly variable
hvg_reachable = scglue.graph.reachable_vertices(guidance, rna.var.query("highly_variable").index)
atac.var["highly_variable"] = [item in hvg_reachable for item in atac.var_names]

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)

# to check the integration consistency
dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, guidance
)
_ = sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")
plt.savefig(f"{datadir}/figures/integration_consistency_{date}.png")

# create the subgraph only of the highly variable features
from itertools import chain
guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()

# store the embeddings 
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

rna.obs["modality"] = "RNA"
atac.obs["modality"] = "ATAC"

# create the combined dataset
combined = ad.concat([rna, atac])

sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)


feature_embeddings = glue.encode_graph(guidance_hvf)
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)

rna.varm["X_glue"] = feature_embeddings.reindex(rna.var_names).to_numpy()
atac.varm["X_glue"] = feature_embeddings.reindex(atac.var_names).to_numpy()

rna.var["name"] = rna.var_names
atac.var["name"] = atac.var_names

genes = rna.var.query("highly_variable").index
peaks = atac.var.query("highly_variable").index

# concatenate the feature indices and embeddings of the modalities
features = pd.Index(np.concatenate([rna.var_names, atac.var_names]))
feature_embeddings = np.concatenate([rna.varm["X_glue"], atac.varm["X_glue"]])

# extract skeleton graph, on which the regulatory interference is conducted
skeleton = guidance_hvf.edge_subgraph(
    e for e, attr in dict(guidance_hvf.edges).items()
    if attr["type"] == "dist"
).copy()


skeleton.number_of_edges()
# returns a regulatory graph containing regulatory score, p-value and q-value as edge attribute for feature pairs in the skeleton graph 
reginf = scglue.genomics.regulatory_inference(
    features, feature_embeddings,
    skeleton=skeleton, random_state=0,
    alternative="greater" 
)

reginf.number_of_edges()

# returns the biadjacency matrix of the regulatory graph, row order: genes, column order: peaks, weight: edge data key
from networkx.algorithms.bipartite import biadjacency_matrix
res = biadjacency_matrix(reginf, genes, peaks, weight="score", dtype=np.float32)
res2 = biadjacency_matrix(reginf, genes, peaks, weight="pval", dtype=np.float32)
dist = biadjacency_matrix(reginf, genes, peaks, weight="dist", dtype=np.float32)

# save as a sparse coo matrix
window = biadjacency_matrix(
    reginf, genes, peaks, weight=None, dtype=np.float32
).tocoo()
dist = window.multiply(dist.toarray())

# as datafame
df = pd.DataFrame({
    "dist": dist.data.astype(int),
    "glue": res.data,
    "pval": res2.data
})

DIST_BINS = [0, 25, 50, 75, 100, 125, 150]  


#  Some check-ups:
def make_dist_bins(dist, bins):
    r"""
    ``bins`` are in KB
    """
    labels = [f"{bins[i]}-{bins[i+1]} kb" for i in range(len(bins) - 1)]
    bins = np.asarray(bins) * 1e3
    return pd.cut(dist, bins, labels=labels, include_lowest=True)
df["dist_bin"] = make_dist_bins(df["dist"], bins=DIST_BINS)



sns.set_style(style="white")
g=sns.displot(df, x="glue",hue="dist_bin",kind="kde",fill=True, aspect = 1.5)
g.set(xlim=(-1, 1))
plt.savefig(os.path.join(f"{datadir}/figures/lncap_distance_plot_{date}.pdf"))

plt.figure()
sns.histplot(df, x="pval", kde=True)
plt.xlim(0,1)
plt.savefig(os.path.join(f"{datadir}/figures/lncap_pval_histplot_{date}.pdf"))

# get peak to gene links
gene2peak_full = reginf.edge_subgraph(
    e for e, attr in dict(reginf.edges).items()
)

links_full = nx.to_pandas_edgelist(gene2peak_full)
# filter based on distance
df_full = links_full[links_full["dist"] > 0]

links_full.to_csv(f"{datadir}/outputs/lncap_peak2gene_links_nonFilt{date}.tsv", sep ="\t")
df_full.to_csv(f"{datadir}/outputs/lncap_peak2gene_links2_nonFilt{date}.tsv"), sep="\t")

# filter with with qval threshold:
gene2peak = reginf.edge_subgraph(
    e for e, attr in dict(reginf.edges).items()
    if attr["qval"] < th)

# get peak to gene links 
links = nx.to_pandas_edgelist(gene2peak)
df = links[links["dist"] > 0]

links.to_csv(f"{datadir}/outputs/lncap_peak2gene_links_"+str(th)+"_qval_{date}.tsv", sep ="\t")
df.to_csv(os.path.joinf"{datadir}/outputs/lncap_peak2gene_links2_"+str(th)+"_qval_{date}.tsv"), sep="\t")

# convert distance to logarithmic
links["dist_log"] = np.log10(links["dist"] + 1)
#sns.histplot(links, x="dist_log", kde=True)
#sns.histplot(links, x="score", kde=True)

plt.figure()
sns.histplot(links[links["dist"] > 0], x="dist_log", kde=True)
plt.savefig(f"{datadir}/figures/lncap_peak2gene_links_histplot_{date}.pdf")

print("regulatory inference ready")

# save the embeddings, combined object and guidance graphs
rna.write(f"{datadir}/anndata/rna-emb-{date}.h5ad", compression="gzip")
atac.write(f"{datadir}/anndata/atac-emb-{date}.h5ad", compression="gzip")
nx.write_graphml(gene2peak, f"{datadir}/outputs/guidance-"+str(th)+"-qval-{date}.graphml.gz")
combined.write(f"{datadir}/anndata/combined-emb-{date}.h5ad", compression="gzip")
nx.write_graphml(skeleton, f"{datadir}/outputs/guidance-skeleton-{date}.graphml.gz")

