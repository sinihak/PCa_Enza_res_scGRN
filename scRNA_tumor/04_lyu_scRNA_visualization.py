import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
import pickle
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from scipy.stats import spearmanr, zscore, mannwhitneyu
import gseapy as gp
import seaborn as sns
import os
import json
import utils as ut
import random
import h5py
import re
from matplotlib.backends.backend_pdf import PdfPages
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text

from pathlib import Path
datadir = Path.cwd() 



gs = f"{datadir}/data/genesets_to_plot.xlsx"
gene_sets = ut.load_gene_sets_excel(gs)
ar_genes = gene_sets["hallmark_androgen_response"]

gs_to_keep = ["SOX4_regulon", "SOX4_linked_resistance", "ARS_stress_linked_sensitive"]
gene_sets = {key: gene_sets[key] for key in gs_to_keep}


# ignore the ribosomal genes
gene_sets["ARS_stress_linked_sensitive"] =  [gene for gene in gene_sets["ARS_stress_linked_sensitive"] if not re.match(r"^(RPL|RPS)", gene)]
# remove canonical androgen response genes
gene_sets["ARS_stress_linked_sensitive"] = list(set(gene_sets["ARS_stress_linked_sensitive"]) - set(ar_genes))




# remove canonical androgen response genes
# ignore the canonical ar response genes in the sensitive sets
for key, genes in gene_sets.items():
    if "sensitive" in key:
        # remove androgen response genes
        gene_sets[key] = list(set(genes) - set(ar_genes))


def assign_color(row, genes_to_highlight):
    if row["gene"] in genes_to_highlight and row["delta"] > 0:
        return "#D62728"
    elif row["gene"] in genes_to_highlight and row["delta"] < 0:
        return "#1800FF"
    else:
        return "#B0B0B0"


# function to plot the delta rank plot
def plot_rank_delta(
    res,
    genes_to_highlight,
    outfile_name,
    figsize=(8, 6),
    percentile=99
):
    """
    Create a ranked delta plot with highlighted genes of interest.
    Parameters
    ----------
    res : pandas.DataFrame
        Must contain a "delta" column; index should be gene names.
    genes_to_highlight : list or set
        Genes to highlight.
    outfile_name : str
        Filename for the saved PDF.
    figsize : tuple
        Figure size.
    percentile : float
        Percentile for clipping delta values.
    """
    # Prepare dataframe            
    res.columns = ["delta"]
    res["gene"] = res.index # add the gene names from index to a new column
    threshold = np.percentile(np.abs(res["delta"]), percentile) # calculate the plotting thresholds based on given percentile 
    res["delta"] = np.clip(res["delta"], -threshold, threshold) # clip within the tresholds
    res["gene_of_interest"] = res.index.isin(genes_to_highlight)
    print(res["gene_of_interest"])
    res["color"] = res.apply(lambda row: assign_color(row, genes_to_highlight), axis=1)    
    res = res.sort_values("delta", ascending=True).reset_index(drop=True)
    res["rank"] = range(1, len(res) + 1)
    # Plot
    plt.figure(figsize=figsize)
    plt.axhline(0, color="black", linestyle="--", linewidth=0.8, zorder=0)
    # plot the background dots (color grey) first under other dots
    background = res[res["color"] == "#B0B0B0"]
    plt.scatter(background["rank"], background["delta"],   c="#B0B0B0", alpha=0.5, s=30, zorder=1)
    # nex plot the foreground (those dots that are not grey, i.e. our genes of interest)
    foreground = res[res["color"] != "#B0B0B0"]
    plt.scatter(foreground["rank"], foreground["delta"], c=foreground["color"], alpha=0.9, s=50, zorder=2)
    # lastly plot the texts
    texts = [plt.text(r["rank"], r["delta"], r["gene"], fontsize=12, zorder=3)
    for _, r in foreground.iterrows()
    ]
    # adjust texts with arrows so that they don"t overlap
    adjust_text(texts,arrowprops=dict(arrowstyle="-", color="black", lw=1),force_text=4)
    plt.xlabel("Rank")
    plt.ylabel("delta")
    plt.savefig(f"{datadir}/scRNA_tumor/figures/{outfile_name}",dpi=300, bbox_inches="tight", pad_inches=0.5)
    plt.close()


def compute_and_save_mean_delta(data, cells, genes, outpath):
    de = data["de"]
    de_df = pd.DataFrame(de, index=cells, columns=genes)
    mean_delta = de_df.mean(axis=0)
    mean_delta.to_json(outpath, orient="index")
    return mean_delta



data = h5py.File(f"{datadir}/scRNA_tumor/outputs/scRNA_Lyu_et_al_epithelium_lemur_de_HSPC_vs_PC.h5", "r")
cells = pd.read_csv(f"{datadir}/scRNA_tumor/outputs/scRNA_Lyu_cells_lemur.csv", index_col=0).iloc[:, 0]
genes = pd.read_csv(f"{datadir}/scRNA_tumor/outputs/scRNA_Lyu_genes_lemur.csv", index_col=0).iloc[:, 0]
# build delta value dataframe
outpath=f"{datadir}/scRNA_tumor/outputs/scRNA_Lyu_mean_delta_HSPC_vs_PC.json"
mean_delta = compute_and_save_mean_delta(data=data, cells=cells, genes=genes, outpath=outpath)
mean_delta = pd.DataFrame(mean_delta)


# --------------------------
# Lemur rank plot (Fig. 2I) 
# --------------------------

plot_rank_delta(res,genes_to_highlight=sox4_targets, outfile_name="scRNA_Lyu_lemur_rank_plot_HSPC_vs_PC.pdf")


####  Multivariate linear regression analysis

adata_concat = ad.read_h5ad(f"{datadir}/scRNA_tumor/anndata/annotated_adata_obj_refined.h5ad")
# keep the epithelial cells
adata_concat = adata_concat[adata_concat.obs["type"].isin(["PC", "HSPC"])]

adata = adata_concat[adata_concat.obs["cell_type"] == "Epithelial cell"]


auc =  pd.read_csv(f"{datadir}/scRNA_tumor/outputs/auc_scores_Lyu_scRNA_cell_rankings.csv")
auc = auc[auc.index.isin(adata.obs_names)].copy()
adata.obs = pd.concat([adata.obs, auc], axis=1)


sample_df = adata.obs.groupby("sample").agg({
    "SOX4_regulon": "mean",
    "ARS_stress_linked_sensitive":"mean",
    "SOX4_linked_resistance":"mean",
    "Liu_PTEN_loss_UP": "mean",
    "SOX4_linked_resistance_wo_SOX4_regulon":"mean",
    "AR_activity": "mean",
    "type": "first"
}).reset_index()



# -----------------------------------------
# Effect of PTEN and AR activity (Fig. 2J) 
# -----------------------------------------

scores_to_test = ["SOX4_regulon", "ARS_stress_linked_sensitive","SOX4_linked_resistance"]

pdf_path = f"{datadir}/scRNA_tumor/figures/scRNA_Lyu_HSPC_PC_ARS_stress_SOX4_linear_regression.pdf"
with PdfPages(pdf_path) as pdf:
    for score in scores_to_test:
        model = smf.ols(
            f"{score} ~ Liu_PTEN_loss_UP + AR_activity + C(type)",
            data=sample_df
        ).fit(cov_type="HC3") # HC3 for small sample size 
        print(score, model.pvalues, model.rsquared_adj)
        rsquared = model.rsquared_adj
        coef = model.params[1:]  # skip intercept
        conf_int = model.conf_int().iloc[1:]
        print(model.summary())
        plt.figure(figsize=(5, 2.5))
        plt.errorbar(coef, coef.index, xerr=[coef - conf_int[0], conf_int[1] - coef],
                fmt="o", color="#2C1440", capsize=5)
        plt.axvline(0, color="gray", linestyle="--")
        plt.ylabel("Predictor")
        plt.title(f"{score}\nAdj. R² = {rsquared:.2f}") # adjusted R2 penalizes if the some of the variables is irrelevant       
        plt.tight_layout()
        plt.xlim(-4, 4.5)
        pdf.savefig(dpi=300)
        plt.close()


# C(type)[T.HSPC]     0.369095
# Liu_PTEN_loss_UP    0.001270
# AR_activity         0.041502
# dtype: float64 0.6942507441334947


# C(type)[T.HSPC]     0.681006
# Liu_PTEN_loss_UP    0.014372
# AR_activity         0.040077
# dtype: float64 0.7334944559128082

# C(type)[T.HSPC]     7.039994e-02
# Liu_PTEN_loss_UP    5.354703e-03
# AR_activity         9.559637e-04
# dtype: float64 0.5674942966045741


# ----------------------------------------------------------------------------------------------
# Effect of SOX4 regulon, PTEN and AR activity on SOX4-linked resistance (Supplementary Fig. 3H) 
# ----------------------------------------------------------------------------------------------

# Check whether the SOX4-associated resistance depends on the SOX4 regulon 
pdf_path = f"{datadir}/scRNA_tumor/figures/scRNA_Lyu_HSPC_PC_SOX4_linked_resistance_dependency_regression_model.pdf"
with PdfPages(pdf_path) as pdf:
    model = smf.ols(
        "SOX4_linked_resistance_wo_SOX4_regulon ~ SOX4_regulon + Liu_PTEN_loss_UP + AR_activity + C(type)",
        data=sample_df
    ).fit(cov_type="HC3") # HC3 for small sample size 
    print(score, model.pvalues, model.rsquared_adj)
    rsquared = model.rsquared_adj
    coef = model.params[1:]  # skip intercept
    conf_int = model.conf_int().iloc[1:]
    print(model.summary())
    plt.figure(figsize=(5, 2.5))
    plt.errorbar(coef, coef.index, xerr=[coef - conf_int[0], conf_int[1] - coef],
            fmt="o", color="#2C1440", capsize=5)
    plt.axvline(0, color="gray", linestyle="--")
    plt.ylabel("Predictor")
    plt.title(f"{score}\nAdj. R² = {rsquared:.2f}") # adjusted R2 penalizes if the some of the variables is irrelevant       
    plt.tight_layout()
    plt.xlim(-1, 1)
    pdf.savefig(dpi=300)
    plt.close()

# ----------------------------------------------------------------
# SOX4-linked resistance phenotype activity (Supplementary Fig. 3I) 
# -----------------------------------------------------------------

y_options = ["SOX4_linked_resistance"]
groups = list(adata.obs["type"].unique())
stats_list = []
palette={"PC":"#4169E1", "HSPC":"#FF4500"}

pdf_path = f"{datadir}/scRNA_tumor/figures/scRNA_Lyu_HSPC_PC_SOX4_linked_resistance_violin.pdf"
with PdfPages(pdf_path) as pdf:
        for y in y_options:
            plt.figure(figsize=(3, 5))
            sns.violinplot(
                data=adata.obs, 
                y=y, 
                hue="type",
                palette=palette,
                linewidth=1.5,
                split=True,
                inner=None
            )
            plt.xticks(rotation=30)
            plt.ylim(0, 0.15)
            plt.tight_layout()
            pdf.savefig(bbox_inches="tight", dpi=300)
            plt.close()
            vals1 = adata.obs.loc[adata.obs["type"] == groups[0], y].dropna()
            vals2 = adata.obs.loc[adata.obs["type"] == groups[1], y].dropna()
            stat, pval = mannwhitneyu(vals1, vals2, alternative="two-sided")
            pval_decimal = f"{pval:.4f}" 
            stats_list.append({
                "gene set": y,
                "comparison": f"{groups[0]} vs {groups[1]}",
                "p_value": pval,
                "p_value_decimal": pval_decimal 
            })
        stats_df = pd.DataFrame(stats_list)


##### Cell lineage-specific PTEN-association #########

pten_data =  adata.obs.groupby("sample")["Liu_PTEN_loss_UP"].mean().reset_index()

pten_median = pten_data["Liu_PTEN_loss_UP"].median()

# assign the PTEN signature activity level
pten_data["PTEN_signature_activity_level"] = np.where(
    pten_data["Liu_PTEN_loss_UP"] < pten_median, "low", "high"
)

# assign PTEN activity status based on the signature (higher -> PTEN likely lost/low)

# PTEN likely lost
pten_data.loc[
    (pten_data["PTEN_signature_activity_level"] == "high"),
    "inferred_PTEN_activity_status"
] = "low"

pten_data.loc[
    (pten_data["PTEN_signature_activity_level"] == "low"),
    "inferred_PTEN_activity_status"
] = "high"

sample_type = ( 
    adata.obs.groupby("sample")["type"] 
    .first() 
    .reset_index()
    )
    

adata.obs["barcode"] = adata.obs.index
pten_metadata = pd.merge(adata.obs,pten_data,on="sample", how="left", sort=False)
pten_metadata.index = adata.obs.index

adata.obs[["inferred_PTEN_activity_status"]] = pten_metadata[["inferred_PTEN_activity_status"]]


# -------------------------------------------------------------
# SOX4 regulon in epithelial sublineages (Supplementary Fig 3J)
# -------------------------------------------------------------

sox4_regulon = ["SESN3","CD24","KCNMB4","MBNL2","ZMAT1","TUBA1A","MARCKS","ZSWIM5","TMSB4X","C16orf89","ZNF608","MET","SOX4"]
X_sox4 = adata[:, sox4_regulon].X.toarray()
sox4_linkage = sch.linkage(X_sox4.T, method="ward")  # Cluster columns
sox4_dendro = sch.dendrogram(sox4_linkage, no_plot=True)["leaves"]
sox4_order = [sox4_regulon[i] for i in sox4_dendro]


dp4 = sc.pl.DotPlot(adata,var_names=sox4_order, groupby="tumortype_celltype",standard_scale="var",categories_order=order)
dp4.legend(colorbar_title="Mean expression")
dp4.style(cmap="magma", smallest_dot=30, largest_dot=300)
dp4.savefig(f"{datadir}/scRNA_tumor/figures/scRNA_Lyu_tumor_type_refined_epithelial_cell_type_SOX4_regulon_dotplot.pdf",bbox_inches="tight", dpi=300)

# ---------------------------------------
# Cell lineage vs. PTEN effect (Fig. 2K) 
# ---------------------------------------

groups = list(adata.obs["inferred_PTEN_activity_status"].unique())
stats_list = []
palette={"low":"#D2042D", "high":"#0096FF"}
y_options = ["ARS_stress_linked_sensitive","SOX4_linked_resistance"]
hue_order = ["low", "high"]

pdf_path = f"{datadir}/scRNA_tumor/figures/scRNA_Lyu_HSPC_PC_ARS_stress_linked_sensitive_PTEN_status_median_violins.pdf"
with PdfPages(pdf_path) as pdf:
    for pheno in adata.obs["type"].unique():  # loop over unique phenotypes
        phenotype = adata[adata.obs["type"] == pheno]
        for y in y_options:
            plt.figure(figsize=(4, 5))
            # Violinplot (split by AR activity)
            sns.violinplot(
                data=phenotype.obs, 
                x="cell_type_refined", 
                y=y, 
                hue="inferred_PTEN_activity_status",
                hue_order=hue_order,
                linewidth=1.5,
                palette=palette,
                split=True,
                inner=None
            )
            global_median = phenotype.obs[y].median()
            plt.axhline(global_median, color="grey", linestyle="--", linewidth=2)
            plt.xticks(rotation=30)
            plt.title(pheno)
            plt.ylim(0, 0.25)
            plt.tight_layout()
            pdf.savefig(bbox_inches="tight", dpi=300)
            plt.close()
            #  Compute statistics per cell type 
            for cell_type in phenotype.obs["cell_type_refined"].unique():
                subset = phenotype.obs[phenotype.obs["cell_type_refined"] == cell_type]
                vals1 = subset.loc[subset["inferred_PTEN_activity_status"] == groups[0], y].dropna()
                vals2 = subset.loc[subset["inferred_PTEN_activity_status"] == groups[1], y].dropna()
                stat, pval = mannwhitneyu(vals1, vals2, alternative="two-sided")
                pval_decimal = f"{pval:.4f}" 
                stats_list.append({
                    "phenotype": pheno,
                    "cell_type": cell_type,
                    "gene set": y,
                    "comparison": f"{groups[0]} vs {groups[1]}",
                    "p_value": pval,
                    "p_value_decimal": pval_decimal 
                })
            stats_df = pd.DataFrame(stats_list)


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
    sc.tl.rank_genes_groups(rna_temp, groupby=groupby, groups=group, reference=reference, method="wilcoxon", layer="log1p", use_raw=False)
    df = sc.get.rank_genes_groups_df(rna_temp, group=group, pval_cutoff=pval_cutoff)
    df = df.dropna(subset=["logfoldchanges"])
    del rna_temp
    return df


def rank_and_jitter(df):
    """Sort the DataFrame by scores and add jitter."""
    ranked_df = df.sort_values(by="scores", ascending=False)
    ranked_df["scores"] += np.random.normal(0, 1e-5, size=len(ranked_df))
    return ranked_df

def run_gsea(ranked_df, gene_sets="MSigDB_Hallmark_2020"):
    """Run GSEA on a ranked DataFrame with specified gene sets."""
    ranked_list = ranked_df[["names", "scores"]]
    ranked_list.set_index("names", inplace=True)
    gsea_result = gp.prerank(rnk=ranked_list,
                              gene_sets=gene_sets,
                              min_size=5, # keep those pathways/gene sets with at least 5 genes
                              max_size=500, # keep those pathways/gene sets with max 500 genes
                              n_permutations=5000)
    return gsea_result

# ------------------------------------------------------------------------------------------------------------------------
# SOX4-linked resistance vs. ARS- and stress-linked sensitivity in mHSPC club-like and basal cells (Supplementary Fig. 3K) 
# ------------------------------------------------------------------------------------------------------------------------

hspc = adata[adata.obs["type"] == "HSPC"]
hspc_club_basal = hspc[hspc.obs["cell_type_refined"].isin(["basal", "club-like"])]

gs_x = "ARS_stress_linked_sensitive"
gs_y =["SOX4_linked_resistance"]

from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()


#lets visualize the correlation of the two larger gene sets as a scatter
pdf_path = f"{datadir}/scRNA_tumor/figures/scRNA_Lyu_HSPC_ARS_stress_SOX4_club_basal_scatters.pdf"
with PdfPages(pdf_path) as pdf:
    for y in gs_y:
        x_norm = scaler.fit_transform(hspc_club_basal.obs[[gs_x]]).flatten()
        y_norm = scaler.fit_transform(hspc_club_basal.obs[[y]]).flatten()
        fig, ax = plt.subplots(figsize=(6, 5))
        sns.regplot(x=x_norm, y=y_norm, scatter=False, ax=ax, color="#D62728")
        sns.scatterplot(x=x_norm, y=y_norm, ax=ax,linewidth=0.2, s=10, color = "#2C1440")
        # calculate corr. coef and p-value
        corr_coef, p_value = spearmanr(x_norm,y_norm)
        ax.set_title(f"Corr: {corr_coef:.2f}, p: {p_value:.2e}")
        plt.tight_layout()
        pdf.savefig(dpi=300, pad_inches=0.5) 
        plt.close() 



# correlation with the SOX4 regulon
gs_x = "ARS_stress_linked_sensitive"
gs_y =["SOX4_regulon"]
corr_coef, p_value = spearmanr(hspc_club_basal.obs[gs_x],hspc_club_basal.obs[gs_y])
print(corr_coef, p_value)
# 0.5124975730541278 0.0

# ----------------------------------------------------------------------------------------------
#  GSEA of SOX4-linked resistance HIGH vs. LOW club-like and basal cells (Supplementary Table 5) 
# ----------------------------------------------------------------------------------------------

q_low = hspc_club_basal.obs["SOX4_linked_resistance"].quantile(0.25)
q_high = hspc_club_basal.obs["SOX4_linked_resistance"].quantile(0.75)

# Assign groups
def assign_group(x):
    if x <= q_low:
        return "Low_S4LR"
    elif x >= q_high:
        return "High_S4LR"
    else:
        return "Intermediate"

# assign the S4LR group
hspc_club_basal.obs["S4LR_group"] = hspc_club_basal.obs["SOX4_linked_resistance"].copy().apply(assign_group)

degs = compute_DEGs(hspc_club_basal, group=["High_S4LR"], groupby="S4LR_group", reference="Low_S4LR", pval_cutoff=0.01)

# remove ribosomal genes
degs_wo_ribo = degs[~degs["names"].str.match(r"^(RPL|RPS)")] 
degs_ranked = rank_and_jitter(degs_wo_ribo) 

gene_sets = ["GO_Biological_Process_2023", "Reactome_2022"]
gsea_result = run_gsea(degs_ranked, gene_sets=gene_sets)

significant_results = gsea_result.res2d[gsea_result.res2d["FDR q-val"] < 0.05]
# Further filter by number of leading-edge genes >= 5
significant_results = significant_results[
    significant_results["Lead_genes"].str.split(";").apply(len) >= 5
]

significant_results.to_csv(f"{datadir}/scRNA_tumor/outputs/scRNA_Lyu_HSPC_basal_club_high_SOX4_vs_low_GSEA_results_wo_ribosomal_genes.csv")
