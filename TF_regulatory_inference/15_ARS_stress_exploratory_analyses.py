
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import scikit_posthocs as sp
from scipy.stats import zscore, kruskal, mannwhitneyu
import utils as ut
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import linkage, leaves_list
from collections import defaultdict
import random
import os
import re
import pickle
import warnings
warnings.filterwarnings("ignore")

from pathlib import Path
datadir = Path.cwd() 

random.seed(32497)
date = "20240307"

rna = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-w-leiden-phases-tf-regulon-auc-seacells-aliases-{date}.h5ad")
# motif activity matrix:
mAct = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/mAct_matrix_atac_{date}.csv")



# fetch the UPR-associated terms 
enrich_ars_stress = pd.read_csv(f"{datadir}/TF_regulatory_inference/TF_regulatory_inference/outputs/gs_enrichment_stress_ARS.tsv", sep="\t")

upr_terms = [
    'Unfolded Protein Response',
    'Unfolded Protein Response (UPR) R-HSA-381119',
    'ATF4 Activates Genes In Response To Endoplasmic Reticulum Stress R-HSA-380994',
    'PERK Regulates Gene Expression R-HSA-381042'
]

upr_genes = sorted({
    gene
    for _, row in enrich_ars_stress.iterrows()
    if row["Term"] in upr_terms
    for gene in row["Genes"].split(";")
})


# ----------------------------------------
# UPR genes dotplot (Supplementary Fig 2B) 
# ----------------------------------------

# hierarchical clustering of genes
rna_sub = rna[:, upr_genes].copy()
link = linkage(rna_sub.X.T, method='ward')
gene_order = [upr_genes[i] for i in leaves_list(link)]

dp = sc.pl.DotPlot(rna, var_names=gene_order, groupby='dataset', standard_scale='var')
dp.style(cmap='magma', smallest_dot=30, largest_dot=300)
dp.legend(colorbar_title='Mean expression')

dp.savefig(f"{datadir}/TF_regulatory_inference/figures/scRNA_lncap_scRNA_UPR_dotplot.pdf", bbox_inches="tight")
plt.close()


# import gseapy as gp
# reactome_library = gp.get_library(name='Reactome_2022', organism='Human')
# reactome_library['XBP1(S) Activates Chaperone Genes R-HSA-381038']

# ----------------------------------------
# XBP1 targets (Supplementary Fig 2C) 
# ----------------------------------------

auc_data = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/auc_scores_lncap_scRNA_XBP1_targets_20240307.csv")
auc_data["SEACell_alias"] = rna.obs["SEACell_alias"]
auc_cols = ['XBP1s_Activates_Chaperone_Genes_R_HSA_381038']


# Group by SEACell_alias and take mean of the selected regulons
df_means = auc_data.groupby('SEACell_alias')[auc_cols].mean().reset_index()

df_means['dataset'] = df_means['SEACell_alias'].str.split("-").str[0]

# df_means now has one row per SEACell_alias, with the mean of each regulon

numeric_cols = df_means.select_dtypes(include='number').columns
df_means[numeric_cols] = df_means[numeric_cols].apply(lambda x: zscore(x, nan_policy='omit'))

pdf_path = f"{datadir}/TF_regulatory_inference/figures/LNCaP_scRNA_XBPIs_seacell_boxplots.pdf"
with PdfPages(pdf_path) as pdf:
    plt.rcParams.update({'font.size': 5})
    for gene_set in auc_cols:
        plt.figure(figsize=(3, 2))
        # Boxplot
        sns.boxplot(
            data=df_means,
            x='dataset',
            y=gene_set,
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
            y=gene_set,
            color='black',
            size=4,
            alpha=0.8,
            dodge=True   # helps if multiple categories
        )
        plt.ylim(-2, 3.5)
        plt.tight_layout(pad=0.2)
        # Save figure
        pdf.savefig(dpi=300, bbox_inches="tight") 
        plt.close()


results = []


for gene_set in auc_cols:
    # Prepare group data
    datasets = df_means['dataset'].unique()
    groups = [df_means[df_means['dataset'] == ds][gene_set].dropna() for ds in datasets]
    # Kruskal-Wallis test
    kw_stat, kw_p = kruskal(*groups)
    # Post-hoc Mann-Whitney U only if KW is significant
    if kw_p < 0.05:
        dmso_vals = df_means[df_means['dataset'] == 'dmso'][gene_set].dropna()
        for ds in datasets:
            if ds == 'dmso':
                continue
            other_vals = df_means[df_means['dataset'] == ds][gene_set].dropna()
            u_stat, u_p = mannwhitneyu(other_vals, dmso_vals, alternative='two-sided')
            results.append({
                "regulon": gene_set,
                "test": "Mann-Whitney U",
                "comparison": f"{ds} vs dmso",
                "statistic": u_stat,
                "p_value": u_p
            })

# Save results to CSV
results_df = pd.DataFrame(results)
print(results_df)
# >>> print(results_df)
#                                         regulon            test      comparison  statistic   p_value
# 0  XBP1s_Activates_Chaperone_Genes_R_HSA_381038  Mann-Whitney U  enz48h vs dmso       89.0  0.208741
# 1  XBP1s_Activates_Chaperone_Genes_R_HSA_381038  Mann-Whitney U    resA vs dmso       35.0  0.174990
# 2  XBP1s_Activates_Chaperone_Genes_R_HSA_381038  Mann-Whitney U    resB vs dmso       96.0  0.005939

#### TF expression and motif activity analyses ####

# -----------------------------------------------
# Dotplot of TF expression (Supplementary Fig 2D) 
# ----------------------------------------------

plt.rcParams.update({'font.size': 10})

X = rna[:, TFs_to_keep].X
groups = rna.obs["dataset"].values
df = pd.DataFrame(X, columns=TFs_to_keep)
df["group"] = groups

# Mean z-score per group per gene 
group_means = df.groupby("group").mean()

dp1=sc.pl.dotplot(rna, var_names=TFs_to_keep, groupby='dataset', dendrogram=False, show=False, return_fig=True, vcenter=0, vmin=-0.20, vmax=0.20)
dp1.style(cmap='magma', smallest_dot=20, largest_dot=200)
dp1.savefig(f"{datadir}/TF_regulatory_inference/figures/regulon_dotplots/LNCaP_scRNA_dotplot_ARS_stress_TFs_expression.pdf", dpi=300, pad_inches=0.5) 


# process the motif activity data
cols = mAct.columns.to_list()
cols = [col.upper() for col in cols] # had some mouse motif matches for some TFs instead, convert to upper case
mAct.columns = cols

TFs_to_keep = ['FOXO3', 'IRF1', 'SREBF2', 'MYC','SREBF1', 'RREB1', 'ZNF282', 'CREB3', 'ATF4', 'CEBPG', 'HIF1A', 'MSANTD3']
mAct = mAct.set_index('UNNAMED: 0') 
TF_mAct = mAct[TFs_to_keep]
TF_mAct['SEACell_alias'] = atac.obs['SEACell_alias']

mAct_per_seacell = TF_mAct.groupby('SEACell_alias').mean()

# extract dataset info (in this case, this is the dominant sample for the metacell)
mAct_per_seacell['dataset'] = mAct_means.index.str.split("-").str[0]

# ------------------------------------------------------
# ATF4 and CEPGB regulon scatters (Supplementary Fig 2E,F) 
# -------------------------------------------------------

auc_wo_self_loops = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/auc_scores_lncap_scRNA_TF_regulons_wo_self_loops_{date}.csv")
# rna.obs.index.equals(auc_wo_self_loops.index)
# True
auc_wo_self_loops['SEACell_alias'] = rna.obs['SEACell_alias']


genes_of_interest = ['ATF4', 'CEBPG']
regulons_of_interest = ['ATF4_52g', 'CEBPG_37g']


pdf_path = f"{datadir}/TF_regulatory_inference//LNCaP_scRNA_ATF4_CEPBG_expression_regulon_activity_scatter.pdf"
with PdfPages(pdf_path) as pdf:
    for gene, regulon in zip(genes_of_interest, regulons_of_interest):
        # process gex data
        gex = pd.DataFrame(rna[:,gene].X.copy().toarray())
        gex.index = rna.obs.index
        gex.columns = [gene]
        gex['SEACell_alias'] = rna.obs['SEACell_alias']
        gex = gex.groupby('SEACell_alias')[gene].mean().reset_index()
        gex['dataset'] = gex['SEACell_alias'].str.split("-").str[0]
        gex[gene] = (gex[gene] - gex[gene].min()) / (gex[gene].max() - gex[gene].min())
        # process regulon activity
        ract = auc_wo_self_loops[[regulon, 'SEACell_alias']]
        ract = ract.groupby('SEACell_alias')[regulon].mean().reset_index()
        ract[regulon] = (ract[regulon] - ract[regulon].min()) / (ract[regulon].max() - ract[regulon].min())
        df = pd.merge(ract, gex, on = 'SEACell_alias')
        fig, ax = plt.subplots(figsize=(5, 5), dpi=200) # Create a new figure for each scatter plot
        sns.regplot(x=df[gene], y=df[regulon], scatter=False, ax=ax, color = '#5b5b5b')   
        # plot the scatter
        p = sns.scatterplot(x=df[gene], y=df[regulon], 
                        data=df, linewidth=0.2,hue='dataset',palette=ut.dataset_palette,
                        s=100, ax=ax, legend=False, alpha=0.8, zorder=0) # lower z-order, beneath other elements 
        # Calculate correlation
        corr_coef, p_value = spearmanr(df[gene], df[regulon])
        # Set the title with correlation coefficient and p-value
        ax.set_title(f"Corr. Coef.: {corr_coef:.2f}, p-value: {p_value:.2e}")
        ax.set_ylim(-0.1, 1.1)
        ax.set_xlim(-0.1, 1.1)
        # Tighten the layout and save the plot
        plt.tight_layout()
        pdf.savefig(dpi=300)
        plt.close()

# ---------------------------------------------------
# TF motif activity box plots (Supplementary Fig 2G) 
# ---------------------------------------------------

data_types = ['motif_activity']
datas = [mAct_means]
pdf_path = f"{datadir}/TF_regulatory_inference/figures/LNCaP_scRNA_TF_mAct_seacell_boxplots.pdf"
with PdfPages(pdf_path) as pdf:
    for data_type, data in zip(data_types, datas):
        # Z scale:
        numeric_cols = data.select_dtypes(include='number').columns
        data[numeric_cols] = data[numeric_cols].apply(lambda x: zscore(x, nan_policy='omit'))
        for TF in TFs_to_keep:
            plt.figure(figsize=(2, 1.75))
            # Boxplot
            sns.boxplot(
                data=data,
                x='dataset',
                y=TF,
                palette=ut.dataset_palette,  
                width=0.8, 
                showcaps=True,
                fill=False,
                showfliers=False
            )
            sns.swarmplot(
                data=data,
                x='dataset',
                y=TF,
                color='black',
                size=3,
                alpha=0.7,
                dodge=True   
            )
            plt.ylim(-4, 4.1)
            plt.tight_layout()
            plt.title(f'{data_type}')
            # Save figure
            pdf.savefig(dpi=300, bbox_inches="tight") # save the current figure into the PDF
            plt.close()

# ------------------------------------------------------------------------------------------
# MYC, SREBF1, SREBF2 motif activity vs. target gene expression (Supplementary Fig 2H, I, J) 
# ------------------------------------------------------------------------------------------

del mAct_per_seacell['dataset']

def min_max(x):
    return (x - x.min(skipna=True)) / (x.max(skipna=True) - x.min(skipna=True))

# re-define TFs to keep 
TFs_to_keep = ['SREBF2', 'MYC','SREBF1']

# load the regulons to get the target genes
regulons = ut.load_regulons_from_csv(os.path.join(datadir, 'outputs/TF_regulons.csv'))
regulons_sub = {key: regulons[key] for key in TFs_to_keep}
target_genes = list(regulons_sub.values())
# flatten to one list
target_genes = list(set(itertools.chain(*target_genes)))

# shouldn't be strictly necessary, but we'll use the unscaled rna anndata here:
rna_unscaled = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/lncap_scRNA_unscaled_phases_seacells_aliases_{date}.h5ad"))

# Subset the rna data to get only the genes of interest
gene_expression_df = pd.DataFrame(rna_unscaled[:, target_genes].X.toarray(), index=rna_unscaled.obs_names, columns=target_genes)

# Add SEACell information to the gene expression DataFrame
gene_expression_df['SEACell_alias'] = rna_unscaled.obs['SEACell_alias'].values

# Group by SEACell and calculate the mean expression for each gene
mExpr_per_seacell = gene_expression_df.groupby('SEACell_alias').mean()

# min-max to same scale for correlation (again, shouldn't be necessary)
mExpr_per_seacell_minmax = mExpr_per_seacell.apply(min_max)
mAct_per_seacell_minmax = mAct_per_seacell.apply(min_max)


def plot_tf_gene_lollipop(regulons_sub, Expr_df, mAct_df, output_pdf_path):
    mact_df = mAct_df.loc[Expr_df.index]
    with PdfPages(output_pdf_path) as pdf:
        for tf, target_genes in regulons_sub.items():
            tf_values = mact_df[tf]
            correlations, p_values = [], []
            for gene in target_genes:
                gene_values = Expr_df[gene]
                corr, p = spearmanr(tf_values, gene_values)
                correlations.append(corr)
                p_values.append(p)
            df = pd.DataFrame({'gene': target_genes, 'correlation': correlations, 'p_value': p_values})
            df = df.sort_values('correlation', ascending=True)
            fig, ax = plt.subplots(figsize=(max(8, len(df)*0.5), 6))
            ax.stem(df['gene'], df['correlation'], basefmt=" ")
            colors = ['grey' if p >= 0.05 else ('red' if c > 0 else 'blue')
                      for c, p in zip(df['correlation'], df['p_value'])]
            ax.scatter(df['gene'], df['correlation'], color=colors, s=100, zorder=3)
            ax.axhline(0, color='grey', linewidth=1, linestyle='--')
            ax.set_xticklabels(df['gene'], rotation=90)
            ax.set_ylabel('Spearman correlation')
            ax.set_title(f'{tf} correlations')
            ax.set_yticks(np.arange(-1, 1.1, 0.5))
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

# Lollipop plots 
plot_tf_gene_lollipop(regulons_sub, mExpr_per_seacell_minmax, mAct_per_seacell_minmax,
    f"{datadir}/TF_regulatory_inference/figures/TF_motif_activity_expression_lollipop_plots_MYC_SREBF1_SREBF2.pdf")