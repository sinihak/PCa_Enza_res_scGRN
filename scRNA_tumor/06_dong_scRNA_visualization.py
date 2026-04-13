import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import h5py
import random
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import ttest_ind, spearmanr
from matplotlib.backends.backend_pdf import PdfPages
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text
from pathlib import Path
datadir = Path.cwd() 

random.seed(50932)


adata = ad.read_h5ad(f"{datadir}/scRNA_tumor/anndata/scRNA_dong_adata_celltypes_annotated.h5ad")

# read the aucell data
auc =  pd.read_csv(f"{aucdir}/scRNA_tumor/outputs/auc_scores_Dong_scRNA_cell_rankings.csv")
auc = auc[auc.index.isin(adata.obs_names)].copy()

# auc.index.equals(adata.obs.index)
# True
adata.obs = pd.concat([adata.obs, auc], axis=1)

# Select features to visualize 
features=["celltypes","patient", "AR_activity", "NE_identity","SOX4_regulon","SOX4_associated_resistance"]

# Filter to epithelium only
adata_epi = adata[adata.obs["celltypes_broad"] == "epithelium"]


adata_epi.uns["patient_colors"] = plt.cm.tab10_r.colors
adata_epi.uns["sample_colors"] = plt.cm.tab20.colors
adata_epi.uns["celltypes_colors"] =   ["#8C6FA1","#5E81AC","#A3BE8C", "#D08770"]


# -----------------------------------
# UMAPs (Supplementary Fig 4C,D,G,H)
# -----------------------------------

filename = f"{datadir}/figures/Dong_et_al_scRNA_UMAPs.pdf"
from pandas.api.types import is_numeric_dtype
with PdfPages(filename) as pdf:
    for feature in features:
        if is_numeric_dtype(adata_epi.obs[feature]):
            sc.pl.umap(
                adata_epi,
                color=feature,
                show=False,
                cmap="coolwarm",
                vmin="p5",
                vmax="p95",
                size=15
            )
        else:
            sc.pl.umap(
                adata_epi,
                color=feature,
                show=False,
                size=15
            )
        plt.tight_layout()
        pdf.savefig(dpi=300)
        plt.close()


# Annotate NE-like cells based on the NE-identity score
adata.obs["celltypes"] = adata.obs["celltypes"].astype(str)
adata.obs.loc[
    (adata.obs["leiden"].isin(["5", "13"])), "celltypes"
] = "NE_like"

# Rename rest of the tumor/luminal cells as luminal-like
adata.obs["celltypes"] = adata.obs["celltypes"].str.replace("tumor_luminal", "luminal-like")

# --------------------------------------------------------------------------------
# 'SOX4_regulon' and 'SOX4_associated_resistance' violins (Supplementary Fig 4E,F)
# --------------------------------------------------------------------------------

adata.obs['celltypes']
features=['SOX4_regulon','SOX4_associated_resistance']

from matplotlib.pyplot import rc_context

palette =["#8C6FA1","#5E81AC","#A3BE8C", "#D08770"]
with rc_context({"figure.figsize": (4.5, 3)}):
    sc.pl.violin(adata_epi, features, groupby="celltypes",stripplot=False,inner="box", palette=palette)
    plt.savefig(f"{datadir}/scRNA_tumor/figures/Dong_scRNA_celltype_violins.pdf", dpi=300, pad_inches=0.5)


# Calculate proportion of lineages in SOX4 program high and low cells 
df = adata_epi.obs[['SOX4_regulon', 'SOX4_associated_resistance','celltypes']]
# Z scoring
df["SOX4_regulon_Z"] = (df["SOX4_regulon"] - df["SOX4_regulon"].mean()) / df["SOX4_regulon"].std()
df["SOX4_associated_resistance_Z"] = (df["SOX4_associated_resistance"] - df["SOX4_associated_resistance"].mean()) / df["SOX4_associated_resistance"].std()


lineage_cols = ["Basal_Z", "Luminal_Z", "Duct_luminal_Z", "NE_like_Z"]


# ---------------------------------------------------------------------------------------
# Epithelial lineage contribution in SOX4 program high/low cells (Supplementary Fig 4E,F)
# ---------------------------------------------------------------------------------------

pdf_path = f'{datadir}figures/Dong_scRNA_SOX4_cell_lineage_by_high_cells.pdf'
with PdfPages(pdf_path) as pdf:
    for var in group_vars:
        group_name = f"{var}_group"
        # assign group based on z score
        df[group_name] = np.where(df[var] > 0, "High", "Low")
        # Compute proportions per cell type
        prop_df = (
            df
            .groupby(["celltypes", group_name])
            .size()
            .reset_index(name="count")
        )
        prop_df["proportion"] = prop_df.groupby("celltypes")["count"].apply(lambda x: x / x.sum())
        # Contingency table: rows = cell types, columns = High/Low
        cont_table = pd.crosstab(df["celltypes"], df[group_name])
        chi2, p, dof, expected = chi2_contingency(cont_table)
        print(f"\nVariable: {var}")
        print(f"Chi-square statistic: {chi2:.4f}")
        print(f"p-value: {p:.4e}")
        print("Contingency table:")
        print(cont_table)
        # Plot stacked barplot: fraction of High/Low cells per cell type
        plot_df = prop_df.pivot(index="celltypes", columns=group_name, values="proportion")
        plt.figure(figsize=(5,4))
        plot_df.plot(kind="bar", stacked=True, color=["#f05a28", "#3a6da7"])  
        plt.ylabel("Fraction of cells")
        plt.tight_layout()
        pdf.savefig(dpi=300)
        plt.close()


# adata_unscaled.obs.index.equals(adata.obs.index)
# True
adata_unscaled.obs[['leiden','SOX4_associated_resistance', 'SOX4_regulon','celltypes','AR_activity']] = adata.obs[['leiden','SOX4_associated_resistance', 'SOX4_regulon','celltypes','AR_activity']].copy()

# Check correlation with AR activity
adata_unscaled.obs[['AR_activity', 'SOX4_associated_resistance']].corr()
adata_unscaled_lum = adata_unscaled[adata_unscaled.obs['celltypes'] == 'luminal-like']
adata_unscaled_lum.obs[['AR_activity', 'SOX4_associated_resistance']].corr()

# adata_unscaled_lum.obs[['AR_activity', 'SOX4_associated_resistance']].corr()
#                             AR_activity  SOX4_associated_resistance
# AR_activity                     1.00000                     0.25055
# SOX4_associated_resistance      0.25055                     1.00000

# A stricter quantile threshold based on the raw activity score for grouping the luminal cells into high and low groups
q_low = adata_unscaled_lum.obs['SOX4_associated_resistance'].quantile(0.75)
q_high = adata_unscaled_lum.obs['SOX4_associated_resistance'].quantile(0.25)

# Assign groups
def assign_group(x):
    if x <= q_low:
        return 'Low_S4LR'
    elif x >= q_high:
        return 'High_S4LR'
    else:
        return 'Intermediate'

# assign the S4LR group
adata_unscaled_lum.obs['S4LR_group'] = adata_unscaled_lum.obs['SOX4_associated_resistance'].copy().apply(assign_group)

def compute_DEGs(adata, group,groupby,reference, pval_cutoff):
    """
    Computes differential expression genes (DEGs) between specified groups and reference from a given SEACell AnnData object.

    Parameters:
    - adata (AnnData): The AnnData object.
    - group (str): The name of the group to compare against the reference.
    - reference (str): The name of the reference group for comparison.
    - pval_cutoff (float): The p-value cutoff for the DEG results.

    Returns:
    - pd.DataFrame: A DataFrame containing the DEG results.
    """
    rna_temp = adata.copy()
    sc.tl.rank_genes_groups(rna_temp, groupby=groupby, groups=group, reference=reference, method='wilcoxon')
    df = sc.get.rank_genes_groups_df(rna_temp, group=group, pval_cutoff=pval_cutoff)
    df = df.dropna(subset=['logfoldchanges'])
    del rna_temp
    return df


def rank_and_jitter(df):
    """Sort the DataFrame by scores and add jitter."""
    ranked_df = df.sort_values(by='scores', ascending=False)
    ranked_df['scores'] += np.random.normal(0, 1e-5, size=len(ranked_df))
    return ranked_df

def run_gsea(ranked_df, gene_sets='MSigDB_Hallmark_2020'):
    """Run GSEA on a ranked DataFrame with specified gene sets."""
    ranked_list = ranked_df[['names', 'scores']]
    ranked_list.set_index('names', inplace=True)
    gsea_result = gp.prerank(rnk=ranked_list,
                              gene_sets=gene_sets,
                              min_size=5, # keep those pathways/gene sets with at least 5 genes
                              max_size=500, # keep those pathways/gene sets with max 500 genes
                              n_permutations=5000)
    return gsea_result


degs = compute_DEGs(adata_unscaled_lum, group=['High_S4LR'], groupby='S4LR_group', reference='Low_S4LR', pval_cutoff=0.01)

degs_ranked = rank_and_jitter(degs)

# ignore ribosomal genes
degs_ranked = degs_ranked[~degs_ranked['names'].str.contains('^(RPL|RPS)')]

gene_sets = ['MSigDB_Hallmark_2020']
gsea_result = run_gsea(degs_ranked,gene_sets)       
significant_results = gsea_result.res2d[gsea_result.res2d['FDR q-val'] < 0.05]

significant_results.to_csv(f'{datadir}outputs/scRNA_Dong_SOX4_status_high_vs_low_GSEA.csv')

significant_results.loc[
    significant_results['Term'].str.contains('Androgen Response'), 'FDR q-val'
].values[0]


# >>> significant_results
#       Name                                               Term        ES       NES NOM p-val FDR q-val FWER p-val   Tag %  Gene %                                         Lead_genes
# 0  prerank            MSigDB_Hallmark_2020__Androgen Response  0.529818  2.347044       0.0       0.0        0.0   32/67  18.89%  PMEPA1;NKX3-1;ZBTB10;ANKH;ABHD2;TSC22D1;KRT8;Z...
# 1  prerank               MSigDB_Hallmark_2020__Myc Targets V1  0.424325  2.058428       0.0  0.002491      0.005  53/130  19.60%  HSP90AB1;HNRNPU;HSPD1;RACK1;SET;SERBP1;TRIM28;...
# 2  prerank             MSigDB_Hallmark_2020__mTORC1 Signaling  0.406198  1.934435       0.0   0.00631      0.019  41/102  22.81%  SQSTM1;HSP90B1;CALR;HSPD1;PPP1R15A;HSPA5;NUPR1...
# 3  prerank               MSigDB_Hallmark_2020__UV Response Dn  0.442379  1.905094       0.0  0.006227      0.025   40/61  38.96%  CITED2;MAP1B;PLCB4;PTGFR;ZMIZ1;F3;GJA1;PHF3;AK...
# 4  prerank                      MSigDB_Hallmark_2020__Hypoxia  0.412837  1.846459  0.001152  0.009365      0.046   40/87  24.10%  JUN;IRS2;BTG1;PNRC1;PPP1R15A;HSPA5;CITED2;CDKN...
# 5  prerank            MSigDB_Hallmark_2020__Protein Secretion  0.411154  1.755894  0.003681  0.024076      0.137   30/57  31.49%  GNAS;KRT18;TSPAN8;CLTA;KIF1B;RAB2A;CLCN3;RAB14...
# 6  prerank               MSigDB_Hallmark_2020__UV Response Up  0.385619   1.74782  0.001164  0.021918      0.145   31/85  21.38%  DNAJA1;SQSTM1;HNRNPU;FOSB;BTG1;SELENOW;EIF5;RH...
# 7  prerank  MSigDB_Hallmark_2020__TNF-alpha Signaling via ...  0.366277  1.725743       0.0  0.023412      0.173  20/106   6.69%  MARCKS;SQSTM1;PMEPA1;JUN;GADD45B;FOSB;IRS2;BTG...
