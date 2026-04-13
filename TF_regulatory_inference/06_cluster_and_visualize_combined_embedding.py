import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
import os
import pandas as pd
import numpy as np
from matplotlib import rcParams

from itertools import chain
from scipy.io import mmread
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
import re
import random

from pathlib import Path
datadir = Path.cwd() 

dataset_palette = {
    "dmso": "#1D61AE",
    "enz48h": "#E56D04",
    "resB": "#6C086D",
    "resA": "#367250"
}

phase_palette = {
    "G1": "#FF4D4D",
    "S": "#66B3FF",
    "G2M": "#00FF7F"
}

random.seed(36751)

date="20240307"
GRNqval="0.1"
rank="10"

# Do cell cycle scoring
CCgenes = [x.strip() for x in open(f"{datadir}/data/regev_lab_cell_cycle_genes.txt")]
s_genes = CCgenes[:43]
g2m_genes = CCgenes[43:]
CCgenes = [x for x in CCgenes if x in rna.var_names]
sc.tl.score_genes_cell_cycle(rna, s_genes=s_genes, g2m_genes=g2m_genes)

datasets = df["dataset"].unique()


rna = ad.read_h5ad(os.path.join(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-{date}.h5ad"))
combined = ad.read_h5ad(os.path.join(f"{datadir}/TF_regulatory_inference/anndata/combined-emb-{date}.h5ad"))
rnaCC = ad.read_h5ad(os.path.join(f"{datadir}/TF_regulatory_inference/anndata/lncap_rna_CC_phases_{date}.h5ad"))
atac = ad.read_h5ad(os.path.join(f"{datadir}/TF_regulatory_inference/anndata/atac-emb-{date}.h5ad"))

# First, check how many cis-regulatory connections we identified:
peak2gene = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_peak2gene_links2_"+str(GRNqval)+"_qval_{date}.tsv", sep="\t")
# create new, empty column to add the peak and gene info (now "source" and "target" are mixed with peak and gene info)
peak2gene["peak"] = None
peak2gene["gene"] = None

# fetch the peaks and genes from "source" column
for i in range(len(peak2gene)):
    source_value = str(peak2gene.at[i, "source"])
    if source_value.startswith("chr"):
        peak2gene.at[i, "peak"] = source_value
    else:
        peak2gene.at[i, "gene"] = source_value

# fetch the peaks and genes from "target" column
for i in range(len(peak2gene)):
    target_value = str(peak2gene.at[i, "target"])
    if target_value.startswith("chr"):
        peak2gene.at[i, "peak"] = target_value
    else:
        peak2gene.at[i, "gene"] = target_value

peak2gene = peak2gene.drop(columns=["Unnamed: 0", "source", "target"])
# keep unique pairs
peak2gene = peak2gene.drop_duplicates()
print(peak2gene.shape)
# (3916, 10) --> 3916 unique peak - gene pairs

rna.obs["phase"] = rnaCC.obs["phase"]

rna.uns["phase_colors"] =phase_palette

# transfer "phase" from rna to atac using GLUE embedding
scglue.data.transfer_labels(rna, atac, "phase", use_rep="X_glue")

combined.obs["phase"] = None

# Assign "phase" from RNA to "combined" where modality is "RNA"
combined.obs.loc[combined.obs["modality"] == "RNA", "phase"] = rna.obs["phase"]

# Assign "phase" from ATAC to "combined" where modality is "ATAC"
combined.obs.loc[combined.obs["modality"] == "ATAC", "phase"] = atac.obs["phase"]

combined.uns["phase_colors"] = phase_palette
combined.uns["dataset_colors"] = dataset_palette

# # cluster the combined embedding
sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)
sc.tl.leiden(combined)

# source data as alternative (CSV)

# --------------------------
# UMAP of datasets (Fig. 1B) 
# --------------------------

sc.pl.embedding(combined, size=10, color=["dataset"], show=False, palette=dataset_palette,basis="X_umap")
plt.tight_layout()
plt.gcf().savefig(f"{datadir}/TF_regulatory_inference/figures/lncap_combined_object_dataset_umap_{date}.pdf")
plt.close()


# --------------------------------------------
# UMAP of data layers (Supplementary Fig. 1A ) 
# --------------------------------------------

sc.pl.umap(combined, size=10, color=["modality"], show=False)
plt.tight_layout()
plt.gcf().savefig(f"{datadir}/TF_regulatory_inference/figures/lncap_combined_object_modality_umap_{date}.pdf")
plt.close()

modalities = ["RNA", "ATAC"]
# Create a dictionary to store the split AnnData objects
ads = {}

# Iterate over each specified modality and create a subset AnnData
for modality in modalities:
    combined_sub = combined[combined.obs["modality"] == modality].copy()
    ads[modality] = combined_sub

rna_comb = ads["RNA"]
atac_comb = ads["ATAC"]

# ---------------------------------------------------
# UMAP of cell cycle phases (Supplementary Fig. 1A ) 
# ---------------------------------------------------

sc.pl.umap(combined, size=10, color=["phase"], show=False, palette=phase_palette)
plt.tight_layout()
plt.gcf().savefig(f"{datadir}/TF_regulatory_inference/lncap_combined_object_CCphase_umap_{date}.pdf")
plt.close()

rna.obs["leiden"] = rna_comb.obs["leiden"]
rna.uns["leiden"] = rna_comb.uns["leiden"]
rna.obsm["X_umap_comb"] = rna_comb.obsm["X_umap"]

atac.obs["leiden"] = atac_comb.obs["leiden"]
atac.uns["leiden"] = atac_comb.uns["leiden"]
atac.obsm["X_umap_comb"] = atac_comb.obsm["X_umap"]

# ---------------------------------------------------------------------
# Fraction of cells in CC phases across samples (Supplementary Fig. 2A) 
# ---------------------------------------------------------------------

df = pd.DataFrame(rna.obs)
# Group by sample and phase, then count the cells in each group
counts = df.groupby(['dataset', 'phase']).size().reset_index(name='cell_count')
# Pivot the table to have samples as rows, phases as columns, and cell counts as values
pivot_table = counts.pivot_table(index='dataset', columns='phase', values='cell_count', fill_value=0)
# Normalize the counts to get proportions
proportions = pivot_table.div(pivot_table.sum(axis=1), axis=0)
percentages = (proportions * 100).round(2)

df = percentages.reset_index().rename_axis(None, axis=1)

datasets = df['dataset'].unique()

phase_palette = {
    'G1': '#FF4D4D',
    'S': '#00FF7F',
    'G2M': '#66B3FF'
}

pdf_path = os.path.join(datadir, "figures/LNCaP_scRNA_CC_pieplots.pdf")

with PdfPages(pdf_path) as pdf:
    for dataset in df['dataset'].unique():
        sub = df[df['dataset'] == dataset]
        values = sub.drop(columns='dataset').iloc[0]
        colors = [phase_palette[col] for col in values.index]
        plt.figure(figsize=(4,4))
        plt.pie(values, labels=values.index, colors=colors,
                autopct='%1.1f%%', startangle=90)
        plt.title(dataset)
        pdf.savefig()   # save the current figure into the PDF
        plt.close()


## save the outputs

import scipy.io as sio
import scipy.sparse as sp
from scipy.io import mmwrite
mmwrite(f"{datadir}/TF_regulatory_inference/outputs/lncap_scRNA_norm_mtx_" + date + ".mtx", rna_sparse)

sio.mmwrite(f"{datadir}/TF_regulatory_inference/outputs/lncap_scRNA_raw_mtx_" + date + ".mtx"),rna.layers["counts"]
sio.mmwrite(f"{datadir}/outputs/lncap_scRNA_raw_mtx_transposed_" + date + ".mtx",rna.layers["counts"].T)
pd.DataFrame(rna.var_names).to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_scRNA_features_{date}.csv")
pd.DataFrame(rna.obs_names).to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_scRNA_cells_{date}.csv")
pd.DataFrame(rna.obs).to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_scRNA_metadata_{date}.csv")

import scipy.sparse as sp

# Convert the matrix to COO format (sparse)
atac_sparse = sp.coo_matrix(atac.X)

mmwrite(f"{datadir}/TF_regulatory_inference/outputs/lncap_scATAC_norm_mtx_" + date + ".mtx", atac_sparse)
pd.DataFrame(atac.var_names).to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_scATAC_features_{date}.csv")
pd.DataFrame(atac.obs_names).to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_scATAC_cells_{date}.csv")
pd.DataFrame(atac.obs).to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_scATAC_metadata_{date}.csv")

rna_comb.write(f"{datadir}/TF_regulatory_inference/anndata/lncap_scRNA_combined_emb_{date}.h5ad")
del atac_comb.uns["phase_colors"] # this causes an error so let's remove it
atac_comb.write(f"{datadir}/TF_regulatory_inference/anndata/lncap_scATAC_combined_emb_{date}.h5ad")

# X_umap_comb will be used for plotting 
rna.write(f"{datadir}/anndata/TF_regulatory_inference/rna-emb-w-leiden-phases-{date}.h5ad", compression="gzip")
atac.write(f"{datadir}/anndata/TF_regulatory_inference/atac-emb-w-leiden-phases-{date}.h5ad", compression = "gzip")
