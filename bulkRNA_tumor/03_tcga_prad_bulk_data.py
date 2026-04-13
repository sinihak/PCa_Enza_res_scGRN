import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import gseapy as gp
import matplotlib
import pickle
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, mannwhitneyu, spearmanr
import seaborn as sns
import utils as ut
import json
import h5py
import random
from matplotlib.backends.backend_pdf import PdfPages


random.seed(50932)


from pathlib import Path
datadir = Path.cwd()

cn_data = pd.read_csv(f"{datadir}/data/TCGA_PRAD_CopyNumber_Gistic2_all_thresholded_by_genes.tsv", sep = "\t")
clin_data = pd.read_csv(f"{datadir}/data/TCGA_PRAD_sampleMap_FPRAD_clinicalMatrix.tsv", sep="\t")

mat = pd.read_csv(f"{datadir}/data/TCGA_PRAD_sampleMap_HiSeqV2.tsv", sep="\t")
mat.index = mat["sample"]
del mat["sample"] 

# Filter out genes with very low total counts
min_total_count = 10  
filtered_mat = mat[mat.sum(axis=1) >= min_total_count]

tcga_z = zscore(filtered_mat, axis=1)

tcga_z.to_csv(f"{datadir}/bulkRNA_tissue/outputs/tcga_prad_z_score_counts.csv")

print(f"Original genes: {mat.shape[0]}")
# Original genes: 20530
print(f"Filtered genes: {filtered_mat.shape[0]}")
# Filtered genes: 19674


pten_data = cn_data[cn_data["Gene Symbol"] == "PTEN"].set_index("Gene Symbol")
# Transpose and fix column name
pten_data = pten_data.transpose().rename(columns={"PTEN": "PTEN_CN"})
# Add sample ID column
pten_data["SAMPLE_ID"] = pten_data.index


with open(f"{datadir}/data/processed_gene_sets.json", "r") as fp:
    gene_sets = json.load(fp)

#################################
######### GSVA Analysis ######### 
#################################


gsva_res = gp.gsva(data=filtered_mat,gene_sets=gene_sets,outdir=None,min_size=0,max_size=100000)
gsva_res_df = gsva_res.res2d.pivot(index="Name", columns="Term", values="ES").reset_index(names="SAMPLE_ID")
gsva_res_df  = gsva_res_df.reset_index(drop=True)
gsva_res_meta = pd.merge(gsva_res_df,pten_data, on="SAMPLE_ID", how="inner") #keeping only matching samples


unaltered = [0]
loss = [-1,-2]
gsva_res_meta.loc[gsva_res_meta["PTEN_CN"].isin(loss), "PTEN_status"] = "pten_altered"
gsva_res_meta.loc[gsva_res_meta["PTEN_CN"].isin(unaltered), "PTEN_status"] = "pten_unaltered"

gsva_res_meta["AR_activity"] = pd.to_numeric(gsva_res_meta["AR_activity"], errors="coerce")
gsva_res_meta["AR_activity_quartile"] = pd.qcut(gsva_res_meta["AR_activity"], q=4, labels=["Q1", "Q2", "Q3","Q4"])

gsva_res_meta["Liu_PTEN_loss"] = gsva_res_meta["Liu_PTEN_loss_UP"] - abs(gsva_res_meta["Liu_PTEN_loss_DOWN"]) 

gsva_res_meta.to_csv(f"{datadir}/bulkRNA_tissue/outputs/tcga_prad_gsva_metadata.csv")


# ---------------------------------
# AR dependency point plot (Fig 2E) 
# ---------------------------------
yax_vals = ["SOX4_regulon"]
counts = gsva_res_meta.groupby(["AR_activity_quartile", "PTEN_status"]).size().reset_index(name="count")

order=["Q4", "Q3", "Q2", "Q1"]
colors={"pten_altered":"#923D0A","pten_unaltered":"#CC7A30"}
sns.set_style("ticks")


pdf_file = f"{datadir}/bulkRNA_tissue/figures/tcga_pointplots_ar_dependency_of_sox4.pdf"
with PdfPages(pdf_file) as pdf:
    for yax in yax_vals:
        plt.figure(figsize=(8, 6))
        sns.pointplot(
            data=gsva_res_meta,
            x="AR_activity_quartile",
            y=yax,
            order=order,
            # palette=colors,
            errorbar="se",  
            markers="o",
            capsize=0.1,
            err_kws={"linewidth": 1},
        )
        plt.xlabel("AR activity quartile")
        plt.ylabel(yax)
        plt.ylim(-0.23, 0.2)
        plt.tight_layout()
        plt.grid(True, alpha=0.3)
        pdf.savefig()
        plt.close()




results = []

quartiles = ["Q1", "Q2", "Q3", "Q4"]
reference = "Q4"
# statistics 

for g in yax_vals:
    # Get reference group (Q4)
    ref_vals = pd.to_numeric(
        gsva_res_meta[gsva_res_meta["AR_activity_quartile"] == reference][g],
        errors="coerce"
    )
    for q in quartiles[:3]:  # Q1, Q2, Q3
        comp_vals = pd.to_numeric(
            gsva_res_meta[gsva_res_meta["AR_activity_quartile"] == q][g],
            errors="coerce"
        )
        stat, p = mannwhitneyu(ref_vals, comp_vals, alternative="two-sided")
        results.append({
            "gene_set": g,
            "reference": reference,
            "comparison": q,
            "p_value": p
        })


results_df = pd.DataFrame(results)
print(results_df)
# 0                 SOX4_regulon       Q1           0.104679             0.035417  0.121257
# 2                 SOX4_regulon       Q2           0.118046            -0.040279  0.002857
# 4                 SOX4_regulon       Q3           0.051875            -0.100641  0.010744
# 6                 SOX4_regulon       Q4           0.080307            -0.115107  0.009926


# -----------------------------------------------------------
# AR dependency point plot stratified by PTEN status (Fig 2F) 
# -----------------------------------------------------------
yax_vals = ["SOX4_regulon", "ARS_stress_linked_sensitive"]

pdf_file = f"{datadir}/bulkRNA_tissue/figures/tcga_pointplots_ar_pten_dependency_of_sox4.pdf"
with PdfPages(pdf_file) as pdf:
    for yax in yax_vals:
        plt.figure(figsize=(8, 6))
        sns.pointplot(
            data=gsva_res_meta,
            x="AR_activity_quartile",
            y=yax,
            order=order,
            palette=colors,
            hue="PTEN_status",
            errorbar="se",  
            markers="o",
            capsize=0.1,
            err_kws={"linewidth": 1},
        )
        plt.xlabel("AR activity quartile")
        plt.ylabel(yax)
        plt.ylim(-0.23, 0.2)
        plt.tight_layout()
        plt.grid(True, alpha=0.3)
        pdf.savefig()
        plt.close()


# Loop through each quartile
quartiles = gsva_res_meta["AR_activity_quartile"].cat.categories

results = []

from scipy.stats import ttest_ind, mannwhitneyu


# statistics 
quartiles = ["Q1", "Q2", "Q3","Q4"]
for q in quartiles:
    for yax in yax_vals:
        # Subset data for the quartile
        subset = gsva_res_meta[gsva_res_meta["AR_activity_quartile"] == q]
        # Separate PTEN status groups
        altered = subset[subset["PTEN_status"] == "pten_altered"][yax]
        unaltered = subset[subset["PTEN_status"] == "pten_unaltered"][yax]
        altered = pd.to_numeric(altered, errors="coerce")
        unaltered = pd.to_numeric(unaltered, errors="coerce")
        stat, p = mannwhitneyu(altered, unaltered, alternative="two-sided")
        results.append({
            "gene_set": yax,
            "quartile": q,
            "pten_altered_mean": altered.mean(),
            "pten_unaltered_mean": unaltered.mean(),
            "p_value": p
        })

results_df = pd.DataFrame(results)
print(results_df)
#                       gene_set quartile  pten_altered_mean  pten_unaltered_mean   p_value
# 0                 SOX4_regulon       Q1           0.104679             0.035417  0.121257
# 1  ARS_stress_linked_sensitive       Q1          -0.118126            -0.202044  0.000325
# 2                 SOX4_regulon       Q2           0.118046            -0.040279  0.002857
# 3  ARS_stress_linked_sensitive       Q2          -0.020695            -0.040885  0.430639
# 4                 SOX4_regulon       Q3           0.051875            -0.100641  0.010744
# 5  ARS_stress_linked_sensitive       Q3           0.044354             0.066425  0.496123
# 6                 SOX4_regulon       Q4           0.080307            -0.115107  0.009926
# 7  ARS_stress_linked_sensitive       Q4           0.124050             0.154662  0.115591

gsva_res_meta.to_csv(f"{datadir}/bulkRNA_tissue/outputs/tcga_gsva_metadata.csv")



# ----------------------------------------------------
# PI3K/AKT/mTOR activity violin (Supplementary Fig 3D) 
# ----------------------------------------------------

gs_to_plot = ["hallmark_pi3k_akt_mtor_signaling"]
gene_sets_to_plot = {key: gene_sets[key] for key in gs_to_plot}

with PdfPages(f"{datadir}/bulkRNA_tissue/figures/tcga_gsva_violinplots.pdf") as pdf:
    for key in gs_to_plot:
        fig, ax = plt.subplots(figsize=(3, 5), dpi=300)
        sns.violinplot( 
            data=gsva_res_meta,
            y=key,
            hue="PTEN_status",
            inner = None,
            linewidth=1,
            linecolor="#2c1440",
            fill=True,
            split=True,
            palette={"pten_altered": "#D2042D", "pten_unaltered":"#0096FF"},
            zorder=1,
            scale="width",
            )
        gsva_res_meta[key] = pd.to_numeric(gsva_res_meta[key], errors="coerce")
        altered_vals = gsva_res_meta.loc[gsva_res_meta["PTEN_status"]=="pten_altered", key]
        unaltered_vals = gsva_res_meta.loc[gsva_res_meta["PTEN_status"]=="pten_unaltered", key]
        # Compute mann-whitney
        stat, pval = mannwhitneyu(altered_vals, unaltered_vals, alternative="two-sided")
        # Annotate p-value above plot
        y_max = gsva_res_meta[key].max()
        ax.text(0.5, 1.02, f"p = {pval:.4f}", ha="center", va="bottom", fontsize=10, transform=ax.transAxes)
        plt.xticks(rotation=30)
        ax.set_ylim(-1.1, 1.1)  # default y-limits
        plt.tight_layout()
        pdf.savefig(fig, dpi=300, pad_inches=0.5)  # Save current figure to PDF
        plt.close(fig) # close to save memory


