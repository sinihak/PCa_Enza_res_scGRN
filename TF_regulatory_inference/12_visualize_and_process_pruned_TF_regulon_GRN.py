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
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
import re
import random

from pathlib import Path
datadir = Path.cwd() 


random.seed(36751)
date="20240307"
num="250" # approx. number of SEACells
rna = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-w-leiden-phases-{date}.h5ad")
atac = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/atac-emb-w-leiden-phases-{date}.h5ad")

rna_sea = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/lncap_rna_seacells_integrated_num_{num}_{date}.h5ad")
atac_sea = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/lncap_atac_seacells_integrated_num_{num}_{date}.h5ad")

rna.obs["SEACell"] = rna_sea.obs["SEACell"]
atac.obs["SEACell"] = atac_sea.obs["SEACell"]
atac.obs["SEACell_confidence"] = atac_sea.obs["SEACell_confidence"]

auc = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/auc_scores_lncap_scRNA_TF_regulons_{date}.csv", header=0)

# keep only the regulons with at least 5 target genes
toFilt = auc.columns.str.extract(r"(\d+)g").astype(float).squeeze() >= 5
discarded_tfs1 = [col.split("_")[0].replace(".", "-") for col in auc.columns[~toFilt.values]]

auc = auc.loc[:, toFilt.values]

tfs = [col.split("_")[0].replace(".", "-") for col in auc.columns]

normalized_counts = rna[:, tfs].X
num_cells_expressing = np.sum(normalized_counts > 0, axis=0)

# Calculate the fraction of cells expressing each TF gene
total_cells = normalized_counts.shape[0]

# expressed_fraction = num_cells_expressing / total_cells

# Create a dataframe with TFs and their expression fractions
tf_expression_df = pd.DataFrame({
    "TF": tfs,
    "Fraction_expressed": expressed_fraction
})

# keep tfs expressed in at least x % of cells
tfs = tf_expression_df[tf_expression_df["Fraction_expressed"] >= 0.02]["TF"]
discarded_tfs2 = tf_expression_df[tf_expression_df["Fraction_expressed"] < 0.02]["TF"].tolist()


tfs = tfs.tolist()
cols = [col for col in auc.columns if col.split("_")[0] in tfs]
# Filter the auc dataframe to consist only of the selected TFs
auc = auc[cols]


# X_umap_comb will be used for plotting 

# ### add the AUCell scores to the RNA object ####
rna.obs = pd.concat([rna.obs, auc], axis=1)

rna.write(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-w-leiden-phases-tf-regulon-auc-seacells-{date}.h5ad")
atac.write(f"{datadir}/TF_regulatory_inference/anndata/atac-emb-w-leiden-seacells-phases-{date}.h5ad")

alias_mapping = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_SEACell_alias_sample_CC_leiden_fraction_{date}.csv")
alias_mapping = alias_mapping[["SEACell", "alias"]]
rna_obs = rna.obs.reset_index()
atac_obs = atac.obs.reset_index()

# Merge DataFrames
merged_rna = pd.merge(rna_obs, alias_mapping, on="SEACell", how="left")
merged_atac = pd.merge(atac_obs, alias_mapping, on="SEACell", how="left")

rna.obs = merged_rna.set_index("index")
rna.obs = rna.obs.rename(columns={"alias": "SEACell_alias"})

atac.obs = merged_atac.set_index("x")
atac.obs = atac.obs.rename(columns={"alias": "SEACell_alias"})

# rename the indexes to be the same
rna.obs.index.name = "barcode"
atac.obs.index.name = "barcode"

rna.write(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-w-leiden-phases-tf-regulon-auc-seacells-aliases-{date}.h5ad")
atac.write(f"{datadir}/TF_regulatory_inference/anndata/atac-emb-w-leiden-phases-seacells-aliases-{date}.h5ad")

pd.DataFrame(rna.obs).to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_scRNA_metadata_w_leiden_phases_tf_regulon_auc_seacells_{date}.csv")
pd.DataFrame(atac.obs).to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_scATAC_metadata_w_leiden_phases_seacells_{date}.csv")

####### visualize the seacells #### 

# The function is fetched from # https://github.com/dpeerlab/SEACells/tree/main/SEACells
def plot_2D(
    ad,
    key="X_umap",
    colour_metacells=True,
    title="Metacell Assignments",
    save_as=None,
    show=True,
    cmap="Set2",
    figsize=(5, 5),
    SEACell_size=20,
    cell_size=10,
):
    """Plot 2D visualization of metacells using the embedding provided in "key".
    :param ad: annData containing "Metacells" label in .obs
    :param key: (str) 2D embedding of data. Default: "X_umap"
    :param colour_metacells: (bool) whether to colour cells by metacell assignment. Default: True
    :param title: (str) title for figure
    :param save_as: (str or None) file name to which figure is saved
    :param cmap: (str) matplotlib colormap for metacells. Default: "Set2"
    :param figsize: (int,int) tuple of integers representing figure size
    """
    umap = pd.DataFrame(ad.obsm[key]).set_index(ad.obs_names).join(ad.obs["SEACell"])
    umap["SEACell"] = umap["SEACell"].astype("category")
    mcs = umap.groupby("SEACell").mean().reset_index()
    plt.figure(figsize=figsize)
    if colour_metacells:
        sns.scatterplot(
            x=0, y=1, hue="SEACell", data=umap, s=cell_size, cmap=cmap, legend=None
        )
        sns.scatterplot(
            x=0,
            y=1,
            s=SEACell_size,
            hue="SEACell",
            data=mcs,
            cmap=cmap,
            edgecolor="black",
            linewidth=1.25,
            legend=None,
        )
    else:
        sns.scatterplot(
            x=0, y=1, color="grey", data=umap, s=cell_size, cmap=cmap, legend=None
        )
        sns.scatterplot(
            x=0,
            y=1,
            s=SEACell_size,
            color="red",
            data=mcs,
            cmap=cmap,
            edgecolor="black",
            linewidth=1.25,
            legend=None,
        )
    plt.xlabel(f"{key}-0")
    plt.ylabel(f"{key}-1")
    plt.title(title)
    ax = plt.gca()
    ax.set_axis_off()
    if save_as is not None:
        plt.savefig(save_as, dpi=300, transparent=True)
    if show:
        plt.show()
    plt.close()



rna = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-w-leiden-phases-tf-regulon-auc-seacells-aliases-{date}.h5ad")
atac = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/atac-emb-w-leiden-phases-seacells-aliases-{date}.h5ad")


combined = ad.read_h5ad(f"{datadir}/anndata/combined-emb-{date}.h5ad")

combined.obs["SEACell"] = None
# Assign "SEACell" from RNA to "combined" where modality is "RNA"
combined.obs.loc[combined.obs["modality"] == "RNA", "SEACell"] = rna.obs["SEACell"]

# Assign "SEACell" from ATAC to "combined" where modality is "ATAC"
combined.obs.loc[combined.obs["modality"] == "ATAC", "SEACell"] = atac.obs["SEACell"]


# -----------------------
# SEACells UMAP (Fig. 1C) 
# -----------------------

plot_2D(combined, key="X_umap",colour_metacells=True)
plt.savefig(f"{datadir}/TF_regulatory_inference/figures/lncap_combined_object_seacells_umap_{date}.pdf"), dpi=300)
plt.close()

plot_2D(combined, key="X_umap", colour_metacells=True, save_as=f"{datadir}/TF_regulatory_inference/figures/lncap_combined_object_seacells_umap_{date}.pdf", show=False)

# -----------------------------------------------------------------------
# 'Persist' signature from Taavitsainen et al 2021 (Supplementary Fig 1C) 
# -----------------------------------------------------------------------


genesets = f"{datadir}/data/genesets_to_plot.xlsx"
gene_sets = ut.load_gene_sets_excel(genesets)

keys_to_extract = ['Persist']

for geneset_name, gene_list in extracted_sets.items():
    sc.tl.score_genes(rna, gene_list=gene_list, score_name=geneset_name)


sc.pl.embedding(rna, size=15, color=['Persist'], show=False, basis='X_umap_comb')
plt.tight_layout()
plt.gcf().savefig(f"{datadir}/TF_regulatory_inference/figures/lncap_scRNA_persist_umap_"+date+".pdf")
plt.close()

# process the TF regulon activity matrix into  TF regulon x dataset_phase format

def auc_in_group(adata, dataset_column, phase_column):
    import pandas as pd
    import numpy as np
    # Initialize final DataFrame outside the loop
    final_auc_df = pd.DataFrame()
    datasets = list(adata.obs[dataset_column].unique())
    phases = list(adata.obs[phase_column].unique())
    for dataset in datasets:
        for phase in phases: 
            # Subset the data for the current group
            sub = adata[(adata.obs[dataset_column] == dataset) & (adata.obs[phase_column] == phase)].copy()
            sub_auc = sub.obs[auc_cols].T
            # Calculate mean activity
            mean_auc = sub_auc.mean(axis=1)
            # Create DataFrame for this group with the correct column name
            column_name = f"{dataset}_{phase}"
            auc_df = pd.DataFrame(mean_auc, columns=[column_name])
            # Accumulate results
            final_auc_df = pd.concat([final_auc_df, auc_df], axis=1)
        del sub, sub_auc
    # Return final DataFrames
    return final_auc_df

final_auc_df=auc_in_group(rna, 'dataset', 'phase')

# ---------------------------------------
# TF regulon activities heatmap (Fig. 1D) 
# ----------------------------------------

col_order = ['dmso_G1', 'enz48h_G1', 'resA_G1', 'resB_G1', 
            'dmso_S', 'enz48h_S', 'resA_S', 'resB_S',
            'dmso_G2M', 'enz48h_G2M', 'resA_G2M', 'resB_G2M']

final_auc_df =  final_auc_df[col_order]

# Define new tick positions
new_ticks = np.arange(0, 1.01, 0.25)  # Adjusting to range from 0 to 1 in steps of 0.25
colors = ["#1a53ff", "#ffff99", "#cc0000"]
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

g = sns.clustermap(final_auc_df.T,
                   annot=False,
                   cmap=custom_cmap,
                   fmt=".2f",
                   figsize=(14, 14),
                   cbar_kws={"label": "Score"},
                   standard_scale=1,
                   metric="euclidean",
                   method="complete",
                   row_cluster=False,
                   col_cluster=True)

# Set new ticks on the color bar
g.cax.set_yticks(new_ticks)
g.cax.set_yticklabels([f"{t:.2f}" for t in new_ticks])
g.savefig(f"{datadir}/TF_regulatory_inference/figures/TF_regulon_CC_phase_heatmap.pdf", dpi=300, bbox_inches="tight")


# Visualize activities of selected TF regulons on metacell level
auc_cols = ['MYC_66g', 'ATF4_52g', 'SOX4_32g', 'SREBF1_19g', 'SREBF2_20g', 'CEBPG_38g']

auc_data['SEACell_alias'] = rna.obs['SEACell_alias']
# Group by SEACell_alias and take mean of the selected regulons
df_means = auc_data.groupby('SEACell_alias')[auc_cols].mean().reset_index()


df_means['dataset'] = df_means['SEACell_alias'].str.split("-").str[0]

# df_means now has one row per SEACell_alias, with the mean of each regulon

numeric_cols = df_means.select_dtypes(include='number').columns
df_means[numeric_cols] = df_means[numeric_cols].apply(lambda x: zscore(x, nan_policy='omit'))

pdf_path = os.path.join(datadir, "/TF_regulatory_inference/figures/LNCaP_scRNA_TF_regulon_seacell_boxplots.pdf")

# --------------------------------------------
# TF regulon activities in metacells (Fig. 1E) 
# --------------------------------------------

with PdfPages(pdf_path) as pdf:
    for regulon in auc_cols:
        plt.figure(figsize=(3, 2.5))
        # Boxplot
        sns.boxplot(
            data=df_means,
            x='dataset',
            y=regulon,
            palette=ut.dataset_palette,  
            width=0.4, 
            showcaps=True,
            fill=False,
            showfliers=False
        )
        # Swarmplot
        sns.swarmplot(
            data=df_means,
            x='dataset',
            y=regulon,
            color='black',
            size=4,
            alpha=0.8,
            dodge=True   # helps if multiple categories
        )
        plt.ylim(-2, 4)
        plt.tight_layout(pad=0.2)
        # Save figure
        pdf.savefig(dpi=300, bbox_inches="tight") # save the current figure into the PDF
        plt.close()

# -----------------------------------------------
# SOX4 correlation with its target genes (Fig 1G) 
# -----------------------------------------------

regulons = ut.load_regulons_from_csv(f'{datadir}/TF_regulatory_inference/outputs/TF_regulons.csv')
TFs_to_keep = ['SOX4']
regulons_sub = {key: regulons[key] for key in TFs_to_keep}
target_genes = list(regulons_sub.values())
# flatten to one list
target_genes = list(set(itertools.chain(*target_genes)))
target_genes.append('SOX4')

# Subset the rna data to get only the genes of interest
gene_expression_df = pd.DataFrame(rna_unscaled[:, target_genes].X.toarray(), index=rna_unscaled.obs_names, columns=target_genes)

# Add SEACell information to the gene expression DataFrame
gene_expression_df['SEACell_alias'] = rna_unscaled.obs['SEACell_alias'].values

# Group by SEACell and calculate the mean expression for each gene
mExpr_per_seacell = gene_expression_df.groupby('SEACell_alias').mean()

seacell_id = rna_unscaled.obs.drop_duplicates('SEACell_alias').set_index('SEACell_alias')['SEACell_identity']
# Map the phase to the mean expression DataFrame
mExpr_per_seacell['SEACell_identity'] = mExpr_per_seacell.index.map(seacell_id)


tf = 'SOX4'
# Get the target genes of this TF
target_genes = regulons_sub.get(tf, [])

mExpr_per_seacell

p_values = []
correlations = []
comparisons = []
# Loop through all TFs in the regulons dictionary

plt.rcParams.update({'font.size': 20})
import math
valid_genes = [gene for gene in target_genes if gene != 'SOX4']

plt.rcParams.update({'font.size': 25})

if valid_genes:
    SOX4_expression = mExpr_per_seacell['SOX4']
    correlations = []
    p_values = []
    # Compute correlations with SOX4
    for gene in valid_genes:
        if gene == 'SOX4':
            continue
        target_expression = mExpr_per_seacell[gene]
        corr, p = spearmanr(SOX4_expression, target_expression)
        correlations.append(corr)
        p_values.append(p)
    # Create DataFrame for plotting
    adjusted = multipletests(p_values, method='fdr_bh')
    padj = adjusted[1]
    df = pd.DataFrame({
        'gene': [gene for gene in valid_genes if gene != 'SOX4'],
        'correlation': correlations,
        'p_value': p_values,
        'pAdj': padj
    })
    # Sort by correlation
    df = df.sort_values('correlation', ascending=True)
    # Plot lollipop
    fig, ax = plt.subplots(figsize=(max(8, len(df)*0.5), 6))
    markerline, stemlines, baseline = ax.stem(df['gene'], df['correlation'], basefmt=" ")
    # Color points by significance
    colors = ['grey' if p >= 0.05 else 'red' if c > 0 else 'blue'
              for c, p in zip(df['correlation'], df['pAdj'])]
    ax.scatter(df['gene'], df['correlation'], color=colors, s=100, zorder=3)
    ax.axhline(0, color='gray', linewidth=1, linestyle='--')
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['gene'], rotation=90)
    ax.set_ylabel('Spearman correlation')
    ax.set_yticks(np.arange(-1, 1.1, 0.5))
    plt.tight_layout()

# Save figure
output_path = f"{datadir}/TF_regulatory_inference/figures/{tf}_regulon_expression_lollipop_plot.pdf"
plt.savefig(output_path)
plt.close(fig)
