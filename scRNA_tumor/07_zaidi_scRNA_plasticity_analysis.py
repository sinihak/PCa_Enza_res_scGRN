
import anndata as ad
import pandas as pd
import scanpy as sc
import seaborn as sns
import doubletdetection
import numpy as np
import scipy.io as sio
import random
import gseapy as gp
from gseapy import dotplot
from scipy.stats import entropy, spearmanr, gaussian_kde
from scipy.io import mmread
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


from pathlib import Path
datadir = Path.cwd() 

random.seed(129412)

adata = ad.read_h5ad(f"{datadir}/scRNA_tumor/anndata/zaidi_scRNA_anndata_clustered.h5ad")

auc = pd.read_csv(f"{datadir}/scRNA_tumor/outputs/auc_scores_Zaidi_scRNA_cell_rankings.csv")

aucols = list(auc.columns)
auc = auc[auc.index.isin(adata.obs.index)]
adata.obs.index

adata.obs = pd.concat([adata.obs, auc], axis=1)

features = ['subtype','patient']

subtype_palette = {
    "CSPC": "#4169E1",  
    "CRPC": "#e60000",  
    "NEPC": "#ff8000"   
    }


# -----------------------------------
# Subtype and patient UMAPs (Fig. 3E, Supplementary Fig 4D) 
# -----------------------------------

sc.pl.umap(adata, size=15, color='subtype', show=False, palette=subtype_palette)
plt.tight_layout()
plt.savefig(f"{datadir}/figures/zaidi_subtype_umap.pdf", dpi=300,bbox_inches="tight", pad_inches=0.5)

palette = (
    sns.color_palette("Set2", 8) +
    sns.color_palette("Dark2", 8) +
    sns.color_palette("Set1", 5)
)

sc.pl.umap(adata, size=15, color='patient', show=False, palette=palette)
plt.tight_layout()
plt.savefig(f"{datadir}/figures/zaidi_patient_umap.pdf", dpi=300,bbox_inches="tight", pad_inches=0.5)



# -------------------------------------------------
# Lineage trend across SOX4 activity bins (Fig. 3F) 
# -------------------------------------------------


df_all = adata.obs

df_all["SOX4_bin"] = pd.qcut(
    df_all["SOX4_regulon"],
    q=4,
    labels=["Q1", "Q2", "Q3", "Q4"]
)

df_all["SOX4_bin"] = df_all["SOX4_bin"].cat.reorder_categories(
    ["Q1", "Q2", "Q3", "Q4"], ordered=True
)

trend = df_all.groupby("SOX4_bin")[[
    "Luminal", "Basal", "Duct.luminal", 'NE_identity'
]].mean()

plt.figure(figsize=(5,4))
trend.plot(marker="o", color=["#C48D7C", "#6882A2", "#A3B892", "#8B759B"], linewidth=2)
plt.xlabel("SOX4-activity quartile")
plt.ylabel("Mean lineage score")
plt.tight_layout()
plt.savefig(f"{datadir}/figures/zaidi_sox4_lineage_trend_all.pdf",
            dpi=300, bbox_inches="tight", pad_inches=0.5)
plt.close()



# ------------------------------------------------------------------
# Density plot of SOX4 regulon activity across tumor types (Fig. 3F) 
# ------------------------------------------------------------------


df["SOX4_tail"] = df["SOX4_regulon"] >= 0.2

tail_frac = pd.crosstab(
    df["subtype"],
    df["SOX4_tail"],
    normalize="index"
)

plt.figure(figsize=(5,5))
sns.kdeplot(
    data=df,
    x="SOX4_regulon",
    hue="subtype",
    fill=False,
    common_norm=False,
    # alpha=0.6,
    palette=subtype_palette
)
# Add tail cutoff
plt.axvline(0.2, linestyle="--", color="black", label="SOX4 regulon-high cutoff")
# Shade tail region
plt.axvspan(0.2, df["SOX4_regulon"].max(), color="grey", alpha=0.15)
plt.xlabel("SOX4 activity")
plt.ylabel("Density")
fracs = tail_frac[True]
plt.text(0.40, 3, f"CSPC: {fracs['CSPC']:.2f}", color="#4169E1")
plt.text(0.40, 2.5, f"CRPC: {fracs['CRPC']:.2f}", color="#e60000")
plt.text(0.40, 2.0, f"NEPC: {fracs['NEPC']:.2f}", color="#ff8000")
plt.legend()
plt.tight_layout()
plt.savefig(f"{datadir}/figures/zaidi_SOX4_disease_subtype_kde.pdf",
            dpi=300, bbox_inches="tight", pad_inches=0.5)
plt.close()


# compute cell-level Shannon entropy for lineage specific signature scores to proxy lineage plasticity:


adata_copy = adata 
df = adata.obs.copy()


# CRPC/NEPC only

df= df[adata.obs["subtype"].isin(["CRPC","NEPC"])]

# compute Shannon's entropy
lineage = df[["Luminal", "Basal", "Duct.luminal","NE_identity"]]

# L1 normalization per cell (probability normalization)
p_lin = lineage.div(lineage.sum(axis=1), axis=0)
df["lineage_entropy"] = entropy(p_lin.T)


df["lineage_bin"] = pd.qcut(
    df["lineage_entropy"].copy(),
    q=4,
    labels=["Q1", "Q2", "Q3", "Q4"]
)

df["lineage_bin"] = df["lineage_bin"].cat.reorder_categories(
    ["Q1", "Q2", "Q3", "Q4"],
    ordered=True
)


spearmanr(df['SOX4_regulon'], df['AR_activity'])
# SignificanceResult(statistic=-0.4460101428722681, pvalue=0.0)

spearmanr(df['SOX4_regulon'], df['lineage_entropy'])
# SignificanceResult(statistic=0.46351779468923754, pvalue=0.0)


AR_SOX4 = df.groupby("lineage_bin")[[
    "AR_activity",
    "SOX4_regulon",
]].mean()


# ------------------------------------------------------------
# SOX4 against AR activity across entropy quartiles (Fig. 3G)
# ------------------------------------------------------------


plt.figure(figsize=(5,5))

sns.kdeplot(data=df, x="AR_activity", y="SOX4_regulon", hue='lineage_bin', palette='bwr')
plt.tight_layout()
plt.savefig(f"{datadir}/figures/zaidi_SOX4_AR_entropy_kde.pdf",
            dpi=300, bbox_inches="tight")
plt.close()

#### Gene set enrichment analysis: functional context of SOX4 activity

# use the unscaled data
adata_unscaled = ad.read_h5ad(f"{datadir}/anndata/zaidi_scRNA_anndata_unscaled.h5ad")

auc = auc[auc.index.isin(adata_unscaled.obs.index)]
adata_unscaled.obs = pd.concat([adata_unscaled.obs, auc], axis=1)

adata_cr = adata_unscaled[adata_unscaled.obs.subtype.isin(['CRPC'])]

# convert to dense
X = adata_cr[:, adata_cr.var_names].X.toarray()

# compute correlation
corr = []

sox4 = adata_cr.obs["SOX4_regulon"].values

for i, gene in enumerate(adata_cr.var_names):
    c, _ = spearmanr(X[:, i], sox4)
    corr.append(c)


corr_df = pd.DataFrame({
    "gene": adata_cr.var_names,
    "spearman_corr": corr
})

# remove NaNs
corr_df = corr_df.dropna()

ranked_genes = (
    corr_df[["gene", "spearman_corr"]]
    .dropna()
    .sort_values(by="spearman_corr", ascending=False)
    .reset_index(drop=True)
)


gene_sets = ['GO_Biological_Process_2023']

pre_res = gp.prerank(
    rnk=ranked_genes,
    gene_sets=gene_sets,
    outdir=None
)
res_sox4 = pre_res.res2d
res_sox4 = res_sox4[res_sox4['FDR q-val'] < 0.05]
res_sox4 = res_sox4.sort_values(by='NES', ascending=False)
res_sox4.to_csv(f"{datadir}/outputs/zaidi_sox4_corr_gsea_result.csv")



df = adata_cr.obs

# To compare the SOX4-associated gene set activity, we run GSEA also for corr. between gene expression and lineage plasticity
# compute Shannon entropy for CRPC cells

lineage = df[["Luminal", "Basal", "Duct.luminal"]]

# L1 normalization per cell (probability normalization)
p_lin = lineage.div(lineage.sum(axis=1), axis=0)
df["lineage_entropy"] = entropy(p_lin.T)

entropy = df["lineage_entropy"].copy()

corr_entropy = []

for i, gene in enumerate(adata_cr.var_names):
    c, _ = spearmanr(X[:, i], entropy)
    corr_entropy.append(c)


corr_entropy = pd.Series(corr_entropy, index=adata_cr.var_names)


ranked_entropy = corr_entropy.sort_values(ascending=False)
ranked_entropy = ranked_entropy.dropna()

pre_res_entr = gp.prerank(
    rnk=ranked_entropy,
    gene_sets=gene_sets,
    outdir=None
)

res_entr = pre_res_entr.res2d
res_entr= res_entr[res_entr['FDR q-val'] < 0.05]
res_entr = res_entr.sort_values(by='NES', ascending=False)
res_entr.to_csv(f"{datadir}/outputs/zaidi_entropy_corr_gsea_result.csv")

to_plot = [res_entr, res_sox4]


# ---------------------------------------------
# GSEA plots (Fig. 3H and Supplementary Fig 4F)
# ---------------------------------------------

filename = (f"{datadir}/figures/zaidi_gsea_plots_sox4_entropy.pdf")
with PdfPages(filename) as pdf:
    for i in to_plot:
        # Sort and take top 15
        plot_df = (i.sort_values("NES", ascending=False).head(15))
        # Dotplot
        ax = dotplot(plot_df,
             column="FDR q-val",
             cmap='Reds',
             size=5,
             figsize=(4,6), top_term=15)
        fig = ax.figure
        pdf.savefig(dpi=300, bbox_inches="tight")
        plt.close()


