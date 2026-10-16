
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gseapy as gp
import utils as ut
import numpy as np
import random
import json
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy.stats import fisher_exact
from scipy.stats import ttest_ind, spearmanr,mannwhitneyu, kruskal, zscore
from matplotlib.backends.backend_pdf import PdfPages
from statsmodels.stats.outliers_influence import variance_inflation_factor
from scipy.stats import chi2_contingency
import os
import re
import warnings
warnings.filterwarnings("ignore")

random.seed(5325)

from pathlib import Path
datadir = Path.cwd()

su2c = pd.read_csv(f"{datadir}/bulkRNA_tissue/outputs/su2c_fpkm_polya_processed.csv", index_col=0)
su2c_meta = pd.read_csv(f"{datadir}/bulkRNA_tissue/outputs/su2c_fpkm_polya_cna_status_metadata_processed.csv")
su2c_meta.columns
# in PTEN-altered situation, AR score?

with open(f"{datadir}/data/processed_gene_sets.json", "r") as fp:
    gene_sets = json.load(fp)

# log transform for
su2c = np.log2(su2c + 1)


#################################
######### GSVA Analysis ######### 
################################# 

gsva_res = gp.gsva(data=su2c,gene_sets=gene_sets,outdir=None,min_size=0,max_size=100000)
gsva_res_df = gsva_res.res2d.pivot(index="Name", columns="Term", values="ES").reset_index(names="SAMPLE_ID")
gsva_res_df  = gsva_res_df.reset_index(drop=True)

su2c_meta = su2c_meta.reset_index(drop=True)
gsva_res_df["SAMPLE_ID"] = gsva_res_df["SAMPLE_ID"].astype(str)
su2c_meta["SAMPLE_ID"] = su2c_meta["SAMPLE_ID"].astype(str)

# merge the GSVA data with metadata
gsva_res_meta = pd.merge(su2c_meta,gsva_res_df, on="SAMPLE_ID", how="inner")


gsva_res_meta["Liu_PTEN_loss"] = gsva_res_meta["Liu_PTEN_loss_UP"] - gsva_res_meta["Liu_PTEN_loss_DOWN"]

# Assign AR activity level based on the Z score of the GSVA score
gsva_res_meta["AR_SCORE_Z"] = (gsva_res_meta["AR_SCORE"] - gsva_res_meta["AR_SCORE"].mean()) / gsva_res_meta["AR_SCORE"].std()
gsva_res_meta["AR_ACTIVITY_LEVEL"] = np.where(
    gsva_res_meta["AR_SCORE_Z"] <= 0, "ar_low", "ar_high"
)

# Assign PTEN activity level based on the Z score of the GSVA score
gsva_res_meta["Liu_PTEN_loss_Z"] = (gsva_res_meta["Liu_PTEN_loss"] - gsva_res_meta["Liu_PTEN_loss"].mean()) / gsva_res_meta["Liu_PTEN_loss"].std()
gsva_res_meta["PTEN_ACTIVITY_LEVEL"] = np.where(
    gsva_res_meta["Liu_PTEN_loss_Z"] <= 0, "pten_high", "pten_low"
)

# remove the samples with homozygous PTEN loss based on the CNV inference but with high PTEN activity 
# mismatch may be due to bulk-level CNV measurement + cellular heterogeneity (e.g. clonality)
gsva_res_meta= gsva_res_meta[~((gsva_res_meta["PTEN_status"] == -2) & (gsva_res_meta["PTEN_ACTIVITY_LEVEL"] == "pten_high"))]


gsva_res_meta["PTEN_AR_ACTIVITY_COMBINED"] = gsva_res_meta["PTEN_ACTIVITY_LEVEL"] + "_" + gsva_res_meta["AR_ACTIVITY_LEVEL"]

# This will be used to plot the correlation plots in bulk_RNAseq_tissue_correlation_analysis.py (Fig 3C and 3D)
gsva_res_meta.to_csv(f"{datadir}/bulkRNA_tissue/outputs/su2c_gsva_metadata.csv")


# gsva_res_meta = pd.read_csv(f"{datadir}/bulkRNA_tissue/outputs/su2c_gsva_metadata.csv")

# ----------------------------------------------------------
# Visualize gene set activities (Fig 3B, supplementary Fig 4C) 
# ----------------------------------------------------------

# from scipy.stats import shapiro

# stat, p = shapiro(gsva_res_meta[SOX4_regulon"].dropna())
# p 
# 0.0009061383316293359

# define gene sets to plot
gene_sets = ["ARS_stress_linked_sensitive","SOX4_linked_resistance","SOX4_regulon", "hallmark_pi3k_akt_mtor_signaling"]


# Fig 3B, supplementary Fig 4C
plot_configs = [
      {
        "data":gsva_res_meta,
        "x": "PTEN_AR_ACTIVITY_COMBINED",
        "filename": f"{datadir}/bulkRNA_tissue/figures/su2c_gsva_violinplots_ar_pten_status_signature.pdf",
        "object_name":"ar_pten_status_signature",
    }
]

plt.rcParams.update({"font.size": 10})

def plot_violin(df, gene_sets, x_col,filename,ref_group, object_name, order=None):
    groups = [g for g in df[x_col].unique() if pd.notna(g)]
    stats_list = []
    with PdfPages(filename) as pdf:
        for key in gene_sets:
            df[key] = pd.to_numeric(df[key], errors="coerce") # ensuring that the data is numeric
            plt.figure(figsize=(3, 4))           
            sns.violinplot(data=df, x=x_col, y=key, inner=None, linewidth=1, linecolor="#2c1440", color="#2c1440", fill=False, order=order)
            # Boxplot
            sns.boxplot(data=df, x=x_col, y=key, width=0.2, showcaps=False,
                        boxprops={"facecolor":"white", "zorder":1, "edgecolor":"#2c1440","linewidth":1},
                        whiskerprops={"linewidth":1, "color":"#2c1440"},
                        medianprops={"linewidth":3, "color":"red"},
                        showfliers=False, order=order)
            # Stripplot
            global_median = df[key].median()
            plt.xticks(rotation=30)
            plt.ylim(-1.1, 1.1)
            plt.tight_layout()
            pdf.savefig(dpi=300)
            plt.close()
            if len(groups) == 2:
                # Mann-Whitney U test for 2 groups
                # subset by treatment status
                vals1 = df.loc[df[x_col] == groups[0], key].dropna()
                vals2 = df.loc[df[x_col] == groups[1], key].dropna()
                stat, pval = mannwhitneyu(vals1, vals2, alternative="two-sided")
                stats_list.append({"object_name": object_name,"gene set": key, "comparison": f"{groups[0]} vs {groups[1]}", "p_value": pval})
            else:
                # Kruskal-Wallis for multiple groups, then Mann-Whitney for pairwise comparisons against a reference group
                vals_all = [df.loc[df[x_col] == g, key].dropna() for g in groups]
                stat_kw, p_kw = kruskal(*vals_all)
                stats_list.append({"object_name": object_name, "gene set": key, "comparison": "Kruskal-Wallis", "p_value": p_kw})
                # Pairwise Mann-Whitney vs reference group if KW significant
                if ref_group is not None and p_kw < 0.05:
                    for g in groups:
                        if g == ref_group:
                            continue
                        vals_ref = df.loc[df[x_col] == ref_group, key].dropna()
                        vals_g = df.loc[df[x_col] == g, key].dropna()
                        if len(vals_ref) < 1 or len(vals_g) < 1:
                            continue
                        stat_w, p_w = mannwhitneyu(vals_ref, vals_g, alternative="two-sided")
                        stats_list.append({"gene set": key, "comparison": f"{ref_group} vs {g}", "p_value": p_w})
    stats_df = pd.DataFrame(stats_list)
    return stats_df

stats_list_all = []

col_order = ["pten_high_ar_low", "pten_low_ar_high", "pten_high_ar_high", "pten_low_ar_low"]

for config in plot_configs:
    if "signature" in config["object_name"]:
        order = col_order
    else: 
        order = None
    stats_df  = plot_violin(config["data"], gene_sets, config["x"], config["filename"],ref_group ="pten_high_ar_high", object_name=config["object_name"], order=order)
    stats_df["p_adj"] = multipletests(stats_df["p_value"], method="fdr_bh")[1]
    stats_list_all.append(stats_df)


all_stats= pd.concat(stats_list_all, ignore_index=True)

all_stats.to_csv(f"{datadir}/bulkRNA_tissue/outputs/su2c_violin_plots_statistics.csv")

gsva_res_meta = gsva_res_meta.rename(columns={
    "Duct luminal" : "Duct_luminal"
 })

scores_to_test = ["Basal", "Luminal","NE_identity","Duct luminal"]

plt.rcdefaults()

# ----------------------------------------------------------
# Regression Analysis (Fig 3E) 
# ----------------------------------------------------------

scores_to_test = ["Basal", "Luminal","NE_identity","Duct_luminal"]
predictors = [
    "AR_variant_ARV", "AR_activity","Liu_PTEN_loss","hallmark_wnt_beta_catenin_signaling","hallmark_pi3k_akt_mtor_signaling","ARS_stress_linked_sensitive","SOX4_linked_resistance","SOX4_regulon" 
]

pdf_path = f"{datadir}/bulkRNA_tissue/figures/SU2C_SOX4_ARS_linear_regression.pdf"
with PdfPages(pdf_path) as pdf:
    for score in scores_to_test:
        formula = f"{score} ~ " + " + ".join(predictors)
        model = smf.ols(
            formula,
            data=gsva_res_meta
        ).fit(cov_type="HC3") # HC3 for small sample size 
        print(score, model.pvalues, model.rsquared_adj)
        rsquared = model.rsquared_adj
        coef = model.params[1:]  # skip intercept
        conf_int = model.conf_int().iloc[1:]
        print(model.summary())
        plt.figure(figsize=(5, 3))
        plt.errorbar(coef, coef.index, xerr=[coef - conf_int[0], conf_int[1] - coef],
                fmt="o", color="#2C1440", capsize=5)
        plt.axvline(0, color="gray", linestyle="--")
        plt.ylabel("Predictor")
        plt.title(f"{score}\nAdj. R² = {rsquared:.2f}") # adjusted R2 penalizes if the some of the variables are irrelevant       
        plt.tight_layout()
        plt.xlim(-1.5, 1.5)
        pdf.savefig(dpi=300)
        plt.close()


# The regression model remarks about large eigenvalue (which might indicate that there are
# strong multicollinearity problems or that the design matrix is singular)
# Let's test for multicollinearity with variance inflation factors:


X = gsva_res_meta[["Liu_PTEN_loss", "AR_activity","SOX4_regulon", "ARS_stress_linked_sensitive","SOX4_linked_resistance", "hallmark_pi3k_akt_mtor_signaling","hallmark_wnt_beta_catenin_signaling","AR_variant_ARV"]]
X = sm.add_constant(X)  # adds intercept

vif_data = pd.DataFrame()
vif_data["feature"] = X.columns
vif_data["VIF"] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
print(vif_data)

#                                feature       VIF
# 0                                const  1.072240
# 1                        Liu_PTEN_loss  1.259894
# 2                          AR_activity  2.195553
# 3                         SOX4_regulon  1.527075
# 4          ARS_stress_linked_sensitive  2.288211
# 5               SOX4_linked_resistance  1.563697
# 6     hallmark_pi3k_akt_mtor_signaling  1.504184
# 7  hallmark_wnt_beta_catenin_signaling  1.235982
# 8                       AR_variant_ARV  1.389747

# VIFs are < 5 (~1-2), so the correlation is not severe and we don't have to apply corrective measures


# ---------------------------------------
# Lineage Contribution Analysis (Fig 3F) 
# ---------------------------------------

# Z score the GSVA values
gsva_res_meta["SOX4_regulon_Z"] = (gsva_res_meta["SOX4_regulon"] - gsva_res_meta["SOX4_regulon"].mean()) / gsva_res_meta["SOX4_regulon"].std()
gsva_res_meta["SOX4_linked_resistance_Z"] = (gsva_res_meta["SOX4_linked_resistance"] - gsva_res_meta["SOX4_linked_resistance"].mean()) / gsva_res_meta["SOX4_linked_resistance"].std()
gsva_res_meta["Basal_Z"] = (gsva_res_meta["Basal"] - gsva_res_meta["Basal"].mean()) / gsva_res_meta["Basal"].std()
gsva_res_meta["Duct_luminal_Z"] = (gsva_res_meta["Duct_luminal"] - gsva_res_meta["Duct_luminal"].mean()) / gsva_res_meta["Duct_luminal"].std()
gsva_res_meta["Luminal_Z"] = (gsva_res_meta["Luminal"] - gsva_res_meta["Luminal"].mean()) / gsva_res_meta["Luminal"].std()
gsva_res_meta["NE_identity_Z"] = (gsva_res_meta["NE_identity"] - gsva_res_meta["NE_identity"].mean()) / gsva_res_meta["NE_identity"].std()


group_vars = [
    "SOX4_linked_resistance_Z",
    "SOX4_regulon_Z"
]

lineage_cols = ["Basal_Z", "Luminal_Z", "Duct_luminal_Z", "NE_identity_Z"]

# Colors for lineages
colors = {
    "Basal_Z": "#E99575",
    "Luminal_Z": "#DB96C0",
    "Duct_luminal_Z": "#95A3C3",
    "NE_identity_Z": "#72B6A1"
}


pdf_path = f"{datadir}/bulkRNA_tissue/figures/SOX4_cell_lineage_stacked_bars.pdf"
with PdfPages(pdf_path) as pdf:
    for var in group_vars:
        group_name = f"{var}_group"
        # assign high and low groups based on the Z score
        gsva_res_meta[group_name] = np.where(gsva_res_meta[var] > 0, "High", "Low")
        # assign dominant cell lineage
        # max value in lineage cols determines the dominant lineage
        gsva_res_meta["dominant_lineage"] = gsva_res_meta[lineage_cols].idxmax(axis=1) 
        # copute proportions
        prop_df = (
            gsva_res_meta
            .groupby([group_name, "dominant_lineage"])
            .size()
            .reset_index(name="count")
        )  
        # Chi-square test
        prop_df["proportion"] = prop_df.groupby(group_name)["count"].apply(lambda x: x / x.sum())
        # Create a contingency table
        cont_table = pd.crosstab(gsva_res_meta[group_name], gsva_res_meta["dominant_lineage"])
        chi2, p, dof, expected = chi2_contingency(cont_table)
        print(f"Chi-square statistic: {chi2:.4f}")
        print(f"p-value: {p:.4e}")
        print("Contingency table:")
        print(cont_table)
        # plot the stacked barplot
        plt.figure(figsize=(2,3))
        plot_df = prop_df.pivot(index=f"{var}_group", columns="dominant_lineage", values="proportion")
        plot_df.plot(kind="bar",stacked=True,
        color=colors,
        )   
        plt.tight_layout()
        pdf.savefig(dpi=300)


# Chi-square statistic: 41.4195
# p-value: 5.3278e-09
# Contingency table:
# dominant_lineage                Basal_Z  Duct_luminal_Z  Luminal_Z  NE_identity_Z
# SOX4_linked_resistance_Z_group                                                   
# High                                 39              37         28             27
# Low                                   7              25         67             19
# <Figure size 200x300 with 0 Axes>
# <Axes: xlabel='SOX4_linked_resistance_Z_group'>
# Chi-square statistic: 22.4823
# p-value: 5.1767e-05
# Contingency table:
# dominant_lineage      Basal_Z  Duct_luminal_Z  Luminal_Z  NE_identity_Z
# SOX4_regulon_Z_group                                                   
# High                       31              28         34             33
# Low                        15              34         61             13


results = []

for var in group_vars:
    group_name = f"{var}_group"
    cont_table = pd.crosstab(gsva_res_meta[group_name], gsva_res_meta["dominant_lineage"])
    for lineage in cont_table.columns:
        contingency = [
            [cont_table.loc["High", lineage], cont_table.loc["Low", lineage]], # high and low within lineage
            [cont_table.loc["High"].sum() - cont_table.loc["High", lineage], # high but not within lineage
            cont_table.loc["Low"].sum() - cont_table.loc["Low", lineage]] # low but not within lineage
        ]
        oddsratio, pval = fisher_exact(contingency)
        results.append({
            "variable": var,
            "lineage": lineage,
            "odds_ratio": oddsratio,
            "p_value": pval
        })

fisher_results = pd.DataFrame(results)

# adjust p-values for multiple testing
from statsmodels.stats.multitest import multipletests
fisher_results["p_adj"] = multipletests(fisher_results["p_value"], method="fdr_bh")[1]

print(fisher_results)

#                    variable         lineage  odds_ratio       p_value         p_adj
# 0  SOX4_linked_resistance_Z         Basal_Z    6.722050  1.018685e-06  4.074741e-06
# 1  SOX4_linked_resistance_Z  Duct_luminal_Z    1.464255  2.405859e-01  3.207812e-01
# 2  SOX4_linked_resistance_Z       Luminal_Z    0.206927  1.248484e-08  9.987876e-08
# 3  SOX4_linked_resistance_Z   NE_identity_Z    1.352733  4.149340e-01  4.149340e-01
# 4            SOX4_regulon_Z         Basal_Z    2.349474  1.408978e-02  2.254365e-02
# 5            SOX4_regulon_Z  Duct_luminal_Z    0.747899  3.796200e-01  4.149340e-01
# 6            SOX4_regulon_Z       Luminal_Z    0.375624  2.633149e-04  7.021730e-04
# 7            SOX4_regulon_Z   NE_identity_Z    3.002481  1.778792e-03  3.557583e-03