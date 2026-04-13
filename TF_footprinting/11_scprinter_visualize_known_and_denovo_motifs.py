import sys
import os
import random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import random
import matplotlib as mpl
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from matplotlib.backends.backend_pdf import PdfPages
from adjustText import adjust_text
import anndata as ad
import scanpy as sc
import json
import csv
import re
import utils as ut
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
datadir = Path.cwd() 

random.seed(553377)

###### Visualize the TF binding score results ######

seacells = pd.read_csv(f"{datadir}/TF_footprinting/data/lncap_scATAC_seacells_formatted_for_scPrinter.csv", index_col=0)
selected_seacells = list(set(seacells["group"]))

# SEACell sample groups
resA_seacells = [i for i in selected_seacells if "resA" in i]
dmso_seacells = [i for i in selected_seacells if "dmso" in i]
enz48h_seacells = [i for i in selected_seacells if "enz48h" in i]
resB_seacells = [i for i in selected_seacells if "resB" in i]



TFBS_scores = pd.read_csv(f"{datadir}data/TFBS_scores_all_seacells_SOX4_regulatory_region_of_interest_150kb_all_folds.csv")
duplicates_mask = TFBS_scores.duplicated(keep=False)  # keep=False marks all duplicates as True
# View all duplicated rows
duplicated_rows = TFBS_scores[duplicates_mask]


TFBS_scores = TFBS_scores.drop("Unnamed: 0", axis=1)

TFBS_scores_long = TFBS_scores.melt(
    id_vars=["chrom", "start", "end", "TF"],
    var_name="SEACell",
    value_name="score"
)

TFBS_scores_long = TFBS_scores_long.drop_duplicates(
    subset=["chrom", "start", "end", "TF", "SEACell"]
)

# extract treatment labels from column names
TFBS_scores_long["dataset"] = (
    TFBS_scores_long["SEACell"]
    .str.replace("TFBS_", "", regex=False)  # remove prefix
    .str.split("-").str[0]                 # keep only condition part
)


TFBS_scores_long["regulatory_region"] = "SOX4"
TFBS_scores_long.to_csv(f"{datadir}plots/TFBS_scores_long_all_seacells_SOX4_regulatory_region_of_interest_150kb.csv")

pivot = TFBS_scores_long.pivot_table(
    index=["TF", "regulatory_region"],
    columns="dataset",
    values="score"
)


# ----------------------------------
# TF binding score barplots (Fig 4C)
# ----------------------------------

rois = ["SOX4"] # define regulatory regions of interest, here SOX4 only

pdf_path = f"{datadir}/TF_footprinting/plots/TFBS_top_TFs_in_regions_of_interest_per_dataset.pdf"
with PdfPages(pdf_path) as pdf:
    for dataset in datasets:
        dataset_color = ut.dataset_palette[dataset] # use dataset-specific palette
        for roi in rois:
        # subset by dataset
            data = TFBS_scores_long[TFBS_scores_long["dataset"] == dataset]
            # keep only the region of interest
            roi_data = data[data["regulatory_region"] == roi]
            # aggregate: mean score per TF
            tf_means = (
                roi_data.groupby("TF")["score"]
                .mean()
                .reset_index()
                .sort_values(by="score", ascending=False)
            )
            tf_means_top = tf_means.head(5)  # top 5
            # plot barplot
            plt.figure(figsize=(4,3))
            sns.barplot(data=tf_means_top, x="score", y="TF",color=dataset_color)
            plt.title(f"{dataset}: {roi} regulatory region",fontsize=9)
            plt.xlabel("Mean TFBS score")
            plt.ylabel("Transcription Factor")
            plt.tight_layout()
            # save current figure to pdf
            pdf.savefig()
            plt.close()



###### Analyze ChromVAR motif scores for the known and de novo motifs ######

chromvar = ad.read_h5ad(f"{datadir}/TF_footprinting/data/chromvar_cisbp.h5ad")

# define motifs of interest
# Here we select NFI motifs and MAF motifs based on the differential accessibility analysis and TF binding scores
# Add also the "master regulator" TFs within our regulons (those with most targets, excluding the cell cycle regulon TFs)
motifs_of_interest = ["NFIC","NFIB", "NFIX","NFIA","FOXA1","AR","SOX4","MYC", "NR3C2", "NR3C1","CTCF", "RREB1","ATF4", "CEBPG", "SREBF1","MAFK","MAF", "MAFG","MAFB"]


def assign_color(row):
    if row["p_adj"] >= 0.05:
        return "gray"          # non-significant
    elif row["t_stat"] > 0:
        return "red"           # significant and positive
    else:
        return "blue"          # significant and negative


date= "20240307"
rna = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-w-leiden-phases-tf-regulon-auc-seacells-aliases-{date}.h5ad")

# Add SEACell_alias
def min_max(x):
    return (x - x.min(skipna=True)) / (x.max(skipna=True) - x.min(skipna=True))

# Select the top seacells based on the regulon activity in the RNA layer

# focus on the non-proliferating (G1) metacells 
rna_sub = rna[rna.obs["SEACell_identity"] == "non-proliferating"]

# read the atac data to match the SEAcells with the chromvar data
atac = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/atac-emb-w-leiden-phases-seacells-aliases-{date}.h5ad")
# chromvar.obs_names.equals(atac.obs_names)

## process known motifs
chromvar.obs["SEACell_alias"] = atac.obs["SEACell_alias"]
chromvar.obs["SEACell_identity"] = atac.obs["SEACell_identity"]
chromvar_g1 = chromvar[chromvar.obs["SEACell_identity"] == "non-proliferating"]
seacell_labels = chromvar_g1.obs["SEACell_alias"]

# process de novo motifs
chromvar_denovo = ad.read_h5ad(f"{datadir}/TF_footprinting/data/denovo/chromvar_de_novo_count.h5ad")
# chromvar_denovo.obs_names.equals(atac.obs_names)
chromvar_denovo.obs["SEACell_alias"] = atac.obs["SEACell_alias"]
chromvar_denovo.obs["SEACell_identity"] = atac.obs["SEACell_identity"]
chromvar_denovo_g1 = chromvar_denovo[chromvar_denovo.obs["SEACell_identity"] == "non-proliferating"]

# # add some de novo motifs of interest 
denovo= ["count_pos_patterns.pattern_0", "count_pos_patterns.pattern_2", "count_pos_patterns.pattern_45","count_pos_patterns.pattern_13","count_pos_patterns.pattern_30","count_neg_patterns.pattern_7","count_pos_patterns.pattern_9"]
motifs_of_interest = motifs_of_interest + denovo

chromvar_datas = [chromvar_g1, chromvar_denovo_g1]
descs = ["known_motifs", "denovo_motifs"]
regulons_to_plot = ["SOX4_32g"]

# -----------------------------------
# ChromVAR t-stat rank plot (Fig 4D)
# -----------------------------------

with PdfPages(f"{datadir}/TF_footprinting/plots/LNCaP_scRNA_scprinter_chromvar_zscore_t_stat_rank_plots_regulon_high_vs_low.pdf") as pdf:
    for regulon in regulons_to_plot:
        ### Compute mean regulon activity per SEACell ###
        df_means = (rna_sub.obs.groupby("SEACell_alias")[regulon].mean().reset_index())
        df_means[f"{regulon}_minmax"] = min_max(df_means[regulon])
        # Order SEACells by scaled activity
        order = (df_means.groupby("SEACell_alias")[f"{regulon}_minmax"].mean().sort_values(ascending=False).index)
        ### BARPLOT ### 
        fig, ax = plt.subplots(figsize=(20, 10))
        sns.barplot(data=df_means,
            x="SEACell_alias",
            y=f"{regulon}_minmax",
            palette="viridis",
            order=order,
            ax=ax
        )
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_title(regulon)
        fig.tight_layout(pad=0.2)
        pdf.savefig(fig, dpi=300)
        plt.close(fig)
        # Identify top/bottom SEACells 
        q25 = df_means[f"{regulon}_minmax"].quantile(0.25)
        q75 = df_means[f"{regulon}_minmax"].quantile(0.75)
        top_seacells = df_means[df_means[f"{regulon}_minmax"] >= q75]["SEACell_alias"]
        down_seacells = df_means[df_means[f"{regulon}_minmax"] <= q25]["SEACell_alias"]
        top_mask = seacell_labels.isin(top_seacells)
        down_mask = seacell_labels.isin(down_seacells)
        ### T‑TEST + RANK PLOT ###
        for chromvar_data, description in zip(chromvar_datas, descs):
            top_vals = chromvar_data.X[top_mask, :]
            down_vals = chromvar_data.X[down_mask, :]
            t_stats, p_values = ttest_ind(
                top_vals, down_vals,
                axis=0, equal_var=False
            )
            _, p_adj, _, _ = multipletests(p_values, method="fdr_bh")
            df = pd.DataFrame({
                "motif": chromvar_data.var_names,
                "t_stat": t_stats,
                "p_value": p_values,
                "p_adj": p_adj
            }).sort_values("t_stat", ascending=False).reset_index(drop=True)
            # Motifs to annotate (motifs of choice and top/bottom 3)
            motifs_to_plot = (motifs_of_interest + list(df["motif"].head(3))+ list(df["motif"].tail(3)))
            # Prepare for plotting
            res = df.sort_values("t_stat", ascending=True).reset_index(drop=True)
            res["significant"] = res["p_adj"] < 0.05
            res["motif_of_interest"] = res["motif"].isin(motifs_to_plot)
            res["color"] = res.apply(assign_color, axis=1)
            res["rank"] = range(1, len(res) + 1)
            ### SCATTER PLOT ### 
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.scatter(res["rank"], res["t_stat"], c=res["color"])
            ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
            # Label selected motifs
            texts = []
            for _, row in res[res["motif_of_interest"]].iterrows():
                texts.append(
                    ax.text(row["rank"], row["t_stat"], row["motif"], fontsize=15)
                )
            adjust_text(
                texts,
                arrowprops=dict(arrowstyle="-", color="black", lw=1),
                force_text=4
            )
            ax.set_xlabel("Rank")
            ax.set_ylabel("t-statistic")
            ax.set_ylim([-65, 80])  
            ax.set_title(f"{regulon} - {description}")
            fig.tight_layout(pad=0.2)
            pdf.savefig(fig, dpi=300)
            plt.close(fig)


# the peaks of interest:
poi_SOX4 = pd.read_csv(f"{datadir}/TF_footprinting/SOX4_regulatory_region_of_interest_150kb_peaks.bed", sep="\t")

# de novo hits
hits_denovo = pd.read_csv(f"{datadir}/TF_footprinting/data/denovo/de_novo_hits_count.tsv", sep="\t")


# generate the interval of the de novo hit locations
hits_denovo[["chr_str", "start_str", "end_str"]] = hits_denovo[["chr", "start", "end"]].astype(str)
hits_denovo["motif_interval"] = hits_denovo[["chr_str", "start_str", "end_str"]].agg("-".join, axis =1)

hits_denovo_sub = hits_denovo[["chr", "start", "end"]]

# rename for pyranges
hits_denovo_sub.columns = ["Chromosome","Start", "End"] 
poi_SOX4.columns = ["Chromosome","Start", "End"] 

# convert to pyranges
gr1 = pr.PyRanges(hits_denovo_sub)
gr2 = pr.PyRanges(poi_SOX4)

# keep overlapping intervals
overlaps = gr1.join(gr2)
overlaps_df = overlaps.df

# generate intervals columns for easier handling
# first for the de novo motifs
overlaps_df[["chr_str", "start_str", "end_str"]] = overlaps_df[["Chromosome", "Start", "End"]].astype(str)
overlaps_df["motif_interval"] = overlaps_df[["chr_str", "start_str", "end_str"]].agg("-".join, axis =1)
# the peaks overlapping the motifs
overlaps_df[["start_b_str", "end_b_str"]] = overlaps_df[["Start_b", "End_b"]].astype(str)
overlaps_df["peak_interval"] = overlaps_df[["chr_str", "start_b_str", "end_b_str"]].agg("-".join, axis =1)

overlaps_df = overlaps_df[["motif_interval", "peak_interval"]]

# de novo NFI and FOXA1 patterns (without the "count" prefix):
moi= ["pos_patterns.pattern_0", "pos_patterns.pattern_2", "pos_patterns.pattern_45", "pos_patterns.pattern_60","pos_patterns.pattern_9", "pos_patterns.pattern_13","pos_patterns.pattern_30","pos_patterns.pattern_51","neg_patterns.pattern_7"]

# keeping the motif intervals overlapping with the SOX4-associated peaks
hits_denovo_SOX4 = hits_denovo[hits_denovo["motif_interval"].isin(overlaps_df["motif_interval"])]
hits_denovo_SOX4_nfi_foxa1 = hits_denovo_SOX4[hits_denovo_SOX4["motif_name"].isin(moi)]

# add the known motif matches
hits_denovo_SOX4_nfi_foxa1.loc[(hits_denovo_SOX4_nfi_foxa1.motif_name=="pos_patterns.pattern_0"),"match0"] = "NFIA"
hits_denovo_SOX4_nfi_foxa1.loc[(hits_denovo_SOX4_nfi_foxa1.motif_name=="pos_patterns.pattern_0"),"match1"] = "NFIC"
hits_denovo_SOX4_nfi_foxa1.loc[(hits_denovo_SOX4_nfi_foxa1.motif_name=="pos_patterns.pattern_0"),"match2"] = "NFIB"

hits_denovo_SOX4_nfi_foxa1.loc[(hits_denovo_SOX4_nfi_foxa1.motif_name=="pos_patterns.pattern_2"),"match0"] = "NFIC"
hits_denovo_SOX4_nfi_foxa1.loc[(hits_denovo_SOX4_nfi_foxa1.motif_name=="pos_patterns.pattern_2"),"match1"] = "NFIB"
hits_denovo_SOX4_nfi_foxa1.loc[(hits_denovo_SOX4_nfi_foxa1.motif_name=="pos_patterns.pattern_2"),"match2"] = "NFIX"

hits_denovo_SOX4_nfi_foxa1.loc[(hits_denovo_SOX4_nfi_foxa1.motif_name=="pos_patterns.pattern_9"),"match0"] = "FOXA2"
hits_denovo_SOX4_nfi_foxa1.loc[(hits_denovo_SOX4_nfi_foxa1.motif_name=="pos_patterns.pattern_9"),"match1"] = "FOXA1"
hits_denovo_SOX4_nfi_foxa1.loc[(hits_denovo_SOX4_nfi_foxa1.motif_name=="pos_patterns.pattern_9"),"match2"] = "FOXB1"

hits_denovo_SOX4_nfi_foxa1.loc[(hits_denovo_SOX4_nfi_foxa1.motif_name=="neg_patterns.pattern_7"),"match0"] = "FOXM1"
hits_denovo_SOX4_nfi_foxa1.loc[(hits_denovo_SOX4_nfi_foxa1.motif_name=="neg_patterns.pattern_7"),"match1"] = "FOXA1"
hits_denovo_SOX4_nfi_foxa1.loc[(hits_denovo_SOX4_nfi_foxa1.motif_name=="neg_patterns.pattern_7"),"match2"] = "FOXK1"

# merge 
hits_denovo_SOX4_nfi_foxa1 = hits_denovo_SOX4_nfi_foxa1.merge(overlaps_df, on="motif_interval", how="left")


### now process the known motif matches
hits_known = pd.read_csv(f"{datadir}/TF_footprinting/data/motif_sites_regulatory_regions_of_interest_150kb_peaks.csv", index_col=0)

# generate the interval of the known motif hit locations
hits_known[["chr_str", "start_str", "end_str"]] = hits_known[["chrom", "start", "end"]].astype(str)
hits_known["motif_interval"] = hits_known[["chr_str", "start_str", "end_str"]].agg("-".join, axis =1)

hits_known_sub = hits_known[["chrom", "start", "end"]]

# rename for pyranges
hits_known_sub.columns = ["Chromosome","Start", "End"] 

gr1 = pr.PyRanges(hits_known_sub)
gr2 = pr.PyRanges(poi_SOX4)

# keep overlapping intervals
overlaps = gr1.join(gr2)
overlaps_df_known = overlaps.df

# motif intervals
overlaps_df_known[["chr_str", "start_str", "end_str"]] = overlaps_df_known[["Chromosome", "Start", "End"]].astype(str)
overlaps_df_known["motif_interval"] = overlaps_df_known[["chr_str", "start_str", "end_str"]].agg("-".join, axis =1)
# the peak intervals overlapping the motifs
overlaps_df_known[["start_b_str", "end_b_str"]] = overlaps_df_known[["Start_b", "End_b"]].astype(str)
overlaps_df_known["peak_interval"] = overlaps_df_known[["chr_str", "start_b_str", "end_b_str"]].agg("-".join, axis =1)

# merge the overlaps and TF information
hits_known_TFs = hits_known[["motif_interval", "TF"]]
overlaps_df_known = overlaps_df_known.merge(hits_known_TFs, on="motif_interval")
TFs_to_keep = ["NFIA", "NFIB", "NFIC","NFIX","FOXA1"]

hits_known_SOX4_nfi_foxa1 = overlaps_df_known[overlaps_df_known["TF"].isin(TFs_to_keep)]

# process known motif hits into a binary peak x TF matrix

known = hits_known_SOX4_nfi_foxa1.copy()
known_binary_matrix = (
    known.assign(value=1) 
    .pivot_table(index="peak_interval", 
    columns="TF", 
    values="value", 
    aggfunc="max", 
    fill_value=0))

# process de novo hits into a binary peak x TF matrix

# keep only subset of columns
hits_denovo_SOX4_nfi_foxa1["matches"] = hits_denovo_SOX4_nfi_foxa1[["match0", "match1", "match2"]].agg("/".join, axis =1)

denovo = hits_denovo_SOX4_nfi_foxa1[["peak_interval","matches"]]


denovo_binary_matrix = (
    denovo.assign(value=1) 
    .pivot_table(index="peak_interval", 
    columns="matches", 
    values="value", 
    aggfunc="max", 
    fill_value=0))


motif_matrix = known_binary_matrix.join(denovo_binary_matrix, how="outer").fillna(0).astype(int)

cols = motif_matrix.columns
# empty dataframe
df = pd.DataFrame(np.zeros((len(cols), len(cols))), index=cols, columns=cols)

# jaccard index
for i in cols:
    for j in cols:
        I = motif_matrix[i]
        print(I)
        J = motif_matrix[j]
        intersection=np.logical_and(I,J)
        union= np.logical_or(I, J)
        df.loc[i, j] = float(intersection.sum()) / float(union.sum())

# --------------------------------------------
# Motif-to-peak jaccard index heatmap (Fig 4E)
# --------------------------------------------

plt.figure(figsize=(5, 5))
ax = sns.heatmap(df, cmap="Reds")
# x_labels
for label in ax.get_yticklabels():
    if label.get_text() in TFs_to_keep:
        label.set_fontweight("bold") # known motifs in bold

# y labels:
for label in ax.get_xticklabels():
    if label.get_text() in TFs_to_keep:
        label.set_fontweight("bold")

        
plt.tight_layout(pad=1)
plt.savefig(f"{datadir}/TF_footprinting/plots/scprinter_FOXA1_NFI_motif_jaccard_heatmap.pdf",bbox_inches="tight")
plt.close()


# The motif match analysis in Supplementary Fig 5A will be done in R with motifmatchR package (scprinter_chromvar_peaks_motif_match_analysis.R)
# Let's visualize NFI expression as UMAPs

# -------------------------------------------------
# NFI gene expression UMAPs (Supplementary Fig 5B)
# -------------------------------------------------


# Define color maps
colors = ["darkblue", "steelblue", "white", "orangered", "darkred"]
cmap = LinearSegmentedColormap.from_list("custom_palette", colors, N=50)

# Function to plot and save UMAP
def plot_umap(data, feature, basis, cmap, filename, vmin=None, vmax=None, data_type=None):
    kwargs = {}
    if vmin is not None:
        kwargs['vmin'] = vmin
    if vmax is not None:
        kwargs['vmax'] = vmax
    ax = sc.pl.embedding(
        data, color=feature, basis=basis, size=20,
        color_map=cmap,
        show=False,
        **kwargs
    )
    # Set title based on data type
    if data_type == 'rna':
        title = f"{feature} expression"
    elif data_type == 'atac':
        title = f"{feature} motif activity"
    else:
        title = feature  
    ax.set_title(title, fontsize=14)
    plt.tight_layout()
    plt.gcf().savefig(filename)
    plt.close()


for feature in TFs_to_keep: 
        expr_path = os.path.join(datadir, f"{datadir}/TF_footprinting//feature_umaps/lncap_scRNA_{feature}_expr_umap_{date}.pdf")
        plot_umap(rna, feature, "X_umap_comb", cmap, expr_path,vmin="p5", vmax="p95", data_type="rna")




# Generate TF data for Fig. 4F

# TFBS data:
TFs = ['NFIA', 'NFIX', 'NFIB', 'NFIC','FOXA1']
scores = TFBS_scores_long[TFBS_scores_long['TF'].isin(TFs)]
scores = scores[['chrom', 'start', 'end', 'TF']]
scores['TF'] = scores['TF'].str.replace(r'^NFI[A-Z]*$', 'NFI', regex=True)


import pyranges as pr

gr_tfbs = pr.PyRanges(scores.rename(columns={
    "chrom":"Chromosome",
    "start":"Start",
    "end":"End"
}))


# AR motif matches:
andr_bed = pd.read_csv(f'{datadir}data/SOX4_region_motif_matches_AR_hocomoco.bed',  sep = '\t', header=None)
andr_bed.columns = ['Chromosome','Start', 'End', 'TF']
andr_bed['TF'] = andr_bed['TF'].str.replace(r'^ANDR.*', 'AR', regex=True)

gr_ar = pr.PyRanges(andr_bed)

## de novo hits
# combine top three significant de novo matches into one column
hits_denovo['known_matches'] = hits_denovo[['match0', 'match1', 'match2']].agg('/'.join, axis=1)

# rename the NFI motif matches as NFI and those that have significant FOX motif matches (incl. FOXA1), rename as FOXA1
def simplify_tf(tf_string):
    parts = tf_string.replace("/", ",").split(",")
    groups = set()
    for p in parts:
        p = p.strip()
        if p.startswith("NFI"):
            groups.add("NFI")
        elif p.startswith("FOX"):
            groups.add("FOXA1")
        else:
            groups.add(p)
    return ",".join(sorted(groups))

denovo["TF"] = denovo["known_matches"].apply(simplify_tf)
denovo = denovo[['chr','start','end','TF']]

denovo_gr = pr.PyRanges(denovo.rename(columns={
    "chr":"Chromosome",
    "start":"Start",
    "end":"End"
}))

# combine AR matches, TFBS, de novo match information into one pyranges object 
gr = pr.concat([gr_ar, gr_tfbs,denovo_gr])

merged = gr.cluster(slack=100)  # merge intervals within 100 bp

df = merged.df

grouped = (
    df.groupby(['Chromosome','Cluster'])
      .agg(
          start=('Start','min'), # smallest end position in cluster
          end=('End','max'), # largest end position in cluster
        TFs=("TF", lambda x: "/".join(sorted(set(x)))) # x is TFs in group, set ensures that the TFs are unique, sort alphabetically and combine into one string
      )
      .reset_index() # turn grouped index back to normal
)

del grouped['Cluster'] 

# This will be used to visualize the track in Fig 4F
grouped.to_csv(f'{datadir}data/scRNA_SOX4_regulatory_region_TFs_clean.bed', sep='\t', header=False, index=False)