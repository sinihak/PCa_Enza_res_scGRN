
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gseapy as gp
import utils as ut
import json
import numpy as np
import random
from scipy.stats import ttest_ind
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import mannwhitneyu
from scipy.stats import spearmanr
import os
import re
import warnings
warnings.filterwarnings("ignore")

random.seed(5325)

from pathlib import Path
datadir = Path.cwd()



dat = pd.read_csv(f"{datadir}/data/GSE126078_norm_counts_TPM_GRCh38.p13_NCBI.tsv",sep="\t",index_col=0) # NCBI-generated TPM matrix from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126078

# Replace Entrez IDs with gene symbols
gene_annot = pd.read_csv(f"{datadir}/data/Human.GRCh38.p13.annot.tsv",sep="\t",low_memory=False)
#>>> gene_annot = pd.read_csv(f"{indir}/Human.GRCh38.p13.annot.tsv",sep="\t")
# <stdin>:1: DtypeWarning: Columns (8,9) have mixed types. Specify dtype option on import or set low_memory=False.

gene_mapping_dict = gene_annot.set_index("GeneID")["Symbol"].to_dict()
dat.index = dat.index.map(gene_mapping_dict)

# Normalize by protein coding
protein_coding_genes = list(set(gene_annot[gene_annot["GeneType"] == "protein-coding"]["Symbol"]))
dat = dat[dat.index.isin(protein_coding_genes)]

meta = pd.read_csv(f"{datadir}/data/GSE126078_series_matrix.txt",sep="\t",skiprows=51,index_col=0)
# Extract relevant columns
meta = meta.T[["!Sample_geo_accession","!Sample_characteristics_ch1","!Sample_source_name_ch1"]]

# Drop unused columns
meta = pd.DataFrame(meta.iloc[:,[0,4,5,7]])
column_names = ["GEO_ID","patient","molecular_phenotype","metastatic_site"]
meta.columns = column_names

# Map the sample names onto columns
geo_to_sample_dict = meta.reset_index().set_index("GEO_ID")["index"].to_dict()
dat.columns = dat.columns.map(geo_to_sample_dict)
meta["molecular_phenotype"].value_counts()

CRPC_only = True
if CRPC_only:
    meta = meta[meta.index.str.contains("CRPC")]
    dat = dat[meta.index.tolist()]
    meta["molecular_phenotype"] = pd.Categorical([s.split(": ")[1] for s in meta["molecular_phenotype"]],categories=["ARpos_NEneg","ARpos_NEpos","ARneg_NEpos","ARlow_NEneg","ARneg_NEneg"])

# Drop the "patient: " prefix in the patient column
meta["patient"] = [p.split(" ")[1] for p in meta["patient"]]

meta_index = meta.index.copy()


# Add the mutation metadata
# From labrecque et al, J Clin Invest 2025, https://doi.org/10.1172/JCI186599 Fig 1D
mutation_metadata = pd.read_excel(f"{datadir}/data/labrecque_2019_mutation_data.xlsx") 

cna_categories = ["amplification", "cn_gain", "cnv", "none", "monoallelic_loss","biallelic_loss","missing"]
for c in mutation_metadata:
    if "_STATUS" in c:
        mutation_metadata[c] = pd.Categorical(mutation_metadata[c],categories = cna_categories)
        mutation_metadata[c] = mutation_metadata[c].cat.remove_unused_categories()

meta = pd.merge(meta,mutation_metadata,left_on="patient",right_on="patient",how="left")
meta[mutation_metadata.columns[1:]] = meta[mutation_metadata.columns[1:]].fillna("missing")

meta.index = meta_index.copy()
tpm = dat.copy()
dat_log = np.log2(tpm + 1)

labrecque_z = zscore(dat_log, axis=1)

labrecque_z.to_csv(f"{datadir}/bulkRNA_tissue/outputs/labrecque_z_score_counts.csv")


#################################
######### GSVA Analysis ######### 
#################################

with open(f"{datadir}/data/processed_gene_sets.json", "r") as fp:
    gene_sets = json.load(fp)


gsva_res = gp.gsva(data=dat_log,gene_sets=gene_sets,outdir=None,min_size=0,max_size=100000)

gsva_res_meta.columns
gsva_res_df = gsva_res.res2d.pivot(index="Name", columns="Term", values="ES").reset_index(names="SAMPLE_ID")
gsva_res_meta = pd.merge(meta,gsva_res_df, right_on="SAMPLE_ID",left_index=True)


gsva_res_meta["Liu_PTEN_loss"] = gsva_res_meta["Liu_PTEN_loss_UP"] - gsva_res_meta["Liu_PTEN_loss_DOWN"]


# keep those with known PTEN status only
gsva_res_meta_sub = gsva_res_meta[gsva_res_meta["PTEN_STATUS"] != "missing"]
list(set(gsva_res_meta_sub["PTEN_STATUS"]))


# assign the simple Pten status 
gsva_res_meta_sub["PTEN_STATUS_SIMPLE"] = np.where(
    gsva_res_meta_sub["PTEN_STATUS"].isin(["biallelic_loss", "monoallelic_loss"]),
    "pten_altered",
    "pten_unaltered"
)


# add a key with empty values for looping
key = "Liu_PTEN_loss"
gene_sets[key] =None


# from scipy.stats import shapiro
# stat, p = shapiro(gsva_res_meta["SOX4_regulon"].dropna())
# p
# 0.008696313016116619

# ------------------------------------------------------------
# Visualize gene set activities (Fig 3A, supplementary Fig 4A) 
# ------------------------------------------------------------

gsva_res_meta_sub["Beltran_NEPC"] = gsva_res_meta_sub["Beltran_NEPC_UP"] - abs(gsva_res_meta_sub["Beltran_NEPC_DOWN"])
gs_to_plot = ["ARS_stress_linked_sensitive","SOX4_associated_resistance","SOX4_regulon", "hallmark_pi3k_akt_mtor_signaling", "AR_activity"]
sns.set_style("ticks")

with PdfPages(f"{datadir}/bulkRNA_tissue/figures/labrecque_gsva_boxplots.pdf") as pdf:
    for key in gs_to_plot:
        fig, ax = plt.subplots(figsize=(3, 3), dpi=300)
        sns.boxplot( 
            data=gsva_res_meta_sub,
            y=key,
            hue="PTEN_STATUS_SIMPLE",
            linewidth=1,
            linecolor="black",
            fill=True,
            palette={"pten_altered": "#D2042D", "pten_unaltered":"#0096FF"},
            zorder=1,
            gap=0.2
            )
        gsva_res_meta_sub[key] = pd.to_numeric(gsva_res_meta_sub[key], errors="coerce")
        altered_vals = gsva_res_meta_sub.loc[gsva_res_meta_sub["PTEN_STATUS_SIMPLE"]=="pten_altered", key]
        unaltered_vals = gsva_res_meta_sub.loc[gsva_res_meta_sub["PTEN_STATUS_SIMPLE"]=="pten_unaltered", key]
        # Compute mann-whitney
        stat, pval = mannwhitneyu(altered_vals, unaltered_vals, alternative="two-sided")
        # Annotate p-value above plot
        y_max = gsva_res_meta_sub[key].max()
        ax.text(0.5, 1.02, f"p = {pval:.4f}", ha="center", va="bottom", fontsize=10, transform=ax.transAxes)
        plt.xticks(rotation=30)
        ax.set_ylim(-1.1, 1.1)  # default y-limits
        plt.tight_layout()
        pdf.savefig(fig, dpi=300, pad_inches=0.5)  # Save current figure to PDF
        plt.close(fig) # close to save memory


# Calculate correlations. We can keep the data with all PTEN statuses
# This will be used to plot the correlation plots in bulk_RNAseq_tissue_correlation_analysis.py (Fig 3C and 3D)

gsva_res_meta.to_csv(f"{datadir}/bulkRNA_tissue/outputs/labrecque_gsva_metadata.csv")


