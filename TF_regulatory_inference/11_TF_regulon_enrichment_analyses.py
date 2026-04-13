
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.io import mmread
from matplotlib.backends.backend_pdf import PdfPages
import re
import gseapy as gp
import utils as ut
import pickle
from pathlib import Path
datadir = Path.cwd() 

# Load target genes
target_genes_df = pd.read_csv(f"{datadir}/outputs/regulons_pruned_grn_{date}.csv")
target_genes_df['genes'] = target_genes_df.iloc[:, 1].apply(lambda x: re.findall(r"'(.*?)'", x))

rna = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-w-leiden-phases-tf-regulon-auc-seacells-aliases-{date}.h5ad")

# Create regulons dictionary
regulons = target_genes_df.groupby('TF')['genes'].apply(list).to_dict()
TFs_to_keep = ["ZNF682","E2F8","E2F7","E2F1","E2F2", "YY1", "NFYB", "ZNF669", "ZNF574", "MSANTD3",
"CEBPG", "ATF4", "CREB3", "ZIC2", "IRF1", "SREBF1", "HIF1A", "SOX4", "FOXO3", "MYC","SREBF2", "ZNF282", "RREB1"]

# Filter regulons to keep only the specified TFs
regulons = {k: regulons[k] for k in TFs_to_keep if k in regulons}
# Extract unique genes
genes = list(set(gene for genes_list in regulons.values() for genes in genes_list for gene in genes))

# define TF regulon groups
stress_responses = ['MSANTD3','CEBPG', 'ATF4', 'CREB3','HIF1A']
resistant=['SOX4','ZIC2','ZNF574']
ARS_associated=['FOXO3', 'IRF1','SREBF2','MYC','SREBF1','ZNF282', 'RREB1']
cell_cycle_related =["ZNF682","E2F8","E2F7","E2F1","E2F2", "YY1", "NFYB", "ZNF669"]

# create a combo of stress and ARS
stress_ARS = list(set(stress_responses + ARS_associated))

# make a dictionary of the TF groups
tf_groups = {
    'stress_responses': stress_responses,
    'resistant': resistant,
    'ARS_associated': ARS_associated,
    'cell_cycle_related':cell_cycle_related,
    'stress_ARS': stress_ARS
}

targets = ut.fetch_target_genes_for_groups(flattened_regulons, tf_groups)
ut.save_regulons_to_csv(targets, os.path.join(datadir, '/TF_regulatory_inference/outputs/TF_regulons_groups.csv'))

gene_sets = ['MSigDB_Hallmark_2020','Reactome_2022']

# perform enrichment analysis for TF regulon groups 
ut.enrichment_analysis(datadir, all_genes, targets,gene_sets=gene_sets, padj_th=0.05)
results = ut.load_enrichment_results(datadir, tf_groups)

def modify_term(term):
    # Replace special characters with "_", convert to uppercase, then remove duplicate underscores
    term = re.sub(r"\W", "_", term).upper()  # Replace non-word characters with "_"
    term = re.sub(r"_+", "_", term)          # Replace multiple underscores with a single "_"
    return term.strip("_")                   # Remove leading/trailing underscores if any


results.to_csv(os.path.join(datadir, f"/TF_regulatory_inference/outputs/gs_enrichment_combined.tsv"), sep="\t", index=False)

results["Gene Count"] = results["Genes"].str.split(";").str.len()
# Then filter the DataFrame for rows where Gene Count is greater than 3
filtered_results = results.loc[results["Gene Count"] >= 3].copy()
# Display the filtered DataFrame

filtered_results.to_csv(os.path.join(datadir, f"/TF_regulatory_inference/outputs/gs_enrichment_combined_filtered.tsv"), sep="\t", index=False)
cols_to_filt = ["Old P-value", "Old adjusted P-value"]
filtered_results = filtered_results.drop(cols_to_filt, axis=1)

# Generate dot plots for Hallmark gene sets
# --------------------------------------------------
# Hallmark enrichment dotplot (Supplementary Fig 1E) 
# --------------------------------------------------

filtered_results = filtered_results[filtered_results['TF Group'] != 'stress_ARS']

dotplot_data = ut.prepare_dotplot_data_with_top_terms(filtered_results, gene_sets, top_term=15, padj_th=0.05)

for gene_set in gene_sets_to_plot:
    ut.plot_dotplot(dotplot_data, gene_set,datadir,top_term=15)

# SEACell clustering
# fetch the AUCell data
auc_cols = [col for col in rna.obs.columns if re.match(r".*_\d+g", col)]

auc = rna.obs[auc_cols].copy()
auc['SEACell_alias'] = rna.obs['SEACell_alias']

# Group by `SEACell` and calculate the mean AUCell values for each term
mean_auc_per_seacell = auc.groupby('SEACell_alias').mean()

# --------------------------------------------
# Metacell cluster plot (Supplementary Fig 1F) 
# --------------------------------------------

# Create a new figure and axes for clustering seacells
fig, ax = plt.subplots(figsize=(12, 3))
# Compute the Euclidean distance matrix based on seacells (rows)
distance_matrix = pdist(mean_auc_per_seacell, metric='euclidean')
# Perform hierarchical clustering (Ward's method minimizes variance within clusters)
linkage_matrix = sch.linkage(distance_matrix, method='ward')
# Generate the dendrogram with customizations for seacell clustering
dendro_result = sch.dendrogram(
    linkage_matrix,
    labels=mean_auc_per_seacell.index,  # Use seacell labels (index of rows in the data)
    orientation='top',
    leaf_rotation=90,  # Rotate the leaf labels
    leaf_font_size=10,  # Font size for leaf labels
    ax=ax  # Specify the axes object for the dendrogram
)
# Manually customize the colors of the branches
# `color_list` contains the colors used for each branch
color_list = []
for i in range(len(dendro_result['color_list'])):
    if i % 3 == 0:  # First third of branches
        color_list.append('teal')
    elif i % 3 == 1:  # Second third of branches
        color_list.append('cornflowerblue')
    else:  # Last third of branches
        color_list.append('peru')
# Apply the custom colors back to the plot
for i, d in enumerate(ax.collections):
    d.set_edgecolor(color_list[i])  # Set the custom edge color for the branches
    d.set_linewidth(1)  # Increase the thickness of the branches

# Set axis labels and tick sizes
ax.set_xlabel("Seacell", fontsize=12)
ax.set_ylabel("Distance", fontsize=12)
ax.tick_params(axis='y', labelsize=12)
ax.tick_params(axis='x', labelsize=12)

plt.savefig(f"{datadir}/TF_regulatory_inference/figures/seacell_clustering_regulon_activity_euclidean_ward.pdf", dpi=300, bbox_inches="tight")
plt.close()


# Extract the order of leaves from the dendrogram
leaf_order = dendro_result['leaves']

# Reorder the `alias_mapping` dataframe to match the dendrogram's seacell order
reordered_cell_cycle_data = alias_mapping.loc[mean_auc_per_seacell.index[leaf_order], ['G1', 'S', 'G2M']]

# Create a separate figure for the barplot of cell cycle proportions
fig2, ax2 = plt.subplots(figsize=(12,2))
# Plot the barplot for the cell cycle phases
# We'll use the same leaf_order for the rows in the bar plot
reordered_cell_cycle_data.plot(kind='bar', stacked=True, color=['#FF0000', '#03FF7F', '#73BAFF'], ax=ax2)
# Customize the barplot
ax2.set_xlabel("Seacell", fontsize=12)
ax2.set_ylabel("Proportion", fontsize=12)
ax2.tick_params(axis='x', rotation=90, labelsize=10)  # Rotate x-axis labels for readability
plt.tight_layout()
plt.savefig(f"{datadir}/TF_regulatory_inference/figures/seacell_clustering_with_cell_cycle_barplot.pdf", dpi=300, bbox_inches="tight")


# Now compute the DEGs for the non-proliferating cells

rna_unscaled = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/lncap_scRNA_unscaled_{date}.h5ad")

rna_unscaled.obs['phase'] = rna.obs['phase'].astype("category")
rna_unscaled.obs['SEACell_alias']= rna.obs['SEACell_alias'].astype("category")


rna.obs['SEACell_identity'] = 'non-proliferating'

proliferating = ['enz48h-9', 'resB-6', 'resA-5', 'resA-6', 'dmso-5', 'enz48h-19','resB-15', 
'resB-16', 'enz48h-5', 'enz48h-8', 'resA-10', 'resA-12','resA-16', 
'resB-14', 'enz48h-4', 'resB-10']

rna.obs.loc[rna.obs['SEACell_alias'].isin(proliferating), 'SEACell_identity'] = 'proliferating'
rna.obs['SEACell_identity'] = rna.obs['SEACell_identity'].astype('category')
rna_scaled.write(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-w-leiden-phases-tf-regulon-auc-seacells-aliases-{date}.h5ad")

rna_unscaled.obs['SEACell_identity'] = rna.obs['SEACell_identity']
rna_unscaled.write(f"{datadir}/TF_regulatory_inference/anndata/lncap_scRNA_unscaled_phases_seacells_aliases_{date}.h5ad"))


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
    sc.tl.rank_genes_groups(rna_temp, groupby=groupby, groups=[group], reference=reference, method='wilcoxon')
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


######## First, calculate DE genes for GRaNPA TF importance analysis ########

pval_cutoff=0.05

df1 = compute_DEGs(rna, group='resA',reference='dmso',pval_cutoff=pval_cutoff)
df2 = compute_DEGs(rna, group='resB',reference='dmso',pval_cutoff=pval_cutoff)
df3 = compute_DEGs(rna, group='enz48h',reference='dmso',pval_cutoff=pval_cutoff)


df1['comparison'] = 'resA_vs_dmso'
df2['comparison'] = 'resB_vs_dmso'
df3['comparison'] = 'enz48h_vs_dmso'


dfs = pd.concat([df1, df2, df3], axis=0)
dfs.to_csv(os.path.join(datadir, 'outputs/lncap_scRNA_dataset_de_genes_vs_dmso.csv'))


######## GSEA analysis of the non-proliferating cells ########

# stricter pval threshold for calculating DE genes for GSEA
pval_cutoff=0.01
log2fc_min=0.5

rna = rna[~rna.obs['SEACell_identity'].isin(['proliferating'])].copy()

gene_sets = ['MSigDB_Hallmark_2020']
datasets = rna.obs['dataset'].unique()


def seacells_compute_gsea(degs, gene_sets='MSigDB_Hallmark_2020', fdr_cutoff=0.05):
    """
    Computes GSEA for each dataset and summarizes the top 10 most altered pathways.
    Parameters:
    - degs (pd.Dataframe): precomputed DEGs.
    - gene_sets (str): The gene set database for GSEA (default: 'MSigDB_Hallmark_2020').
    - fdr_cutoff (float): FDR q-value cutoff for filtering significant GSEA results (default: 0.05).
    Returns:
    - summary_df (pd.DataFrame): A DataFrame summarizing the top 10 most altered pathways across metacells.
    - all_gsea_results (dict): Dictionary with metacell IDs as keys and GSEA results as values.
    """
    # Store all GSEA results
    all_gsea_results = {}
    # Get unique datasets
    # Track top pathways
    pathway_counter = Counter()
    # Loop through each dataset
    for dataset, degs in all_degs.items():
        print(f"Processing dataset: {dataset}")
        if dataset == 'dmso':
            continue    # continue here
        degs_ranked = rank_and_jitter(degs)
        # Run GSEA
        gsea_result = run_gsea(degs_ranked, gene_sets=gene_sets)
        all_gsea_results[dataset] = gsea_result
        # Filter for significant pathways
        print('filtering significant')
        significant_results = gsea_result.res2d[gsea_result.res2d['FDR q-val'] < fdr_cutoff]
        print('fetching significant')
        top_pathways = significant_results['Term'].tolist()  # Extract significant pathways
        pathway_counter.update(top_pathways)
    # Summarize the significantly altered pathways
    # This just calculates the frequency of the significantly altered pathways across datasets
    # you can define, e.g. top 10 most common if desired
    print('summarizing')
    summary_df = pd.DataFrame(pathway_counter.most_common(), columns=['Pathway', 'Frequency'])
    return summary_df, all_gsea_results

degs_dmso = compute_DEGs(rna, group='dmso', groupby='dataset', reference='rest', pval_cutoff=pval_cutoff)
degs_resA = compute_DEGs(rna, group='resA', groupby='dataset', reference='dmso', pval_cutoff=pval_cutoff)
degs_resB = compute_DEGs(rna, group='resB', groupby='dataset', reference='dmso', pval_cutoff=pval_cutoff)
degs_enz48h = compute_DEGs(rna, group='enz48h', groupby='dataset', reference='dmso', pval_cutoff=pval_cutoff)

all_degs = {}
all_degs['dmso'] = degs_dmso
all_degs['enz48h'] = degs_enz48h
all_degs['resA'] = degs_resA
all_degs['resB'] = degs_resB


# Save DEGs
with open(f"{datadir}/TF_regulatory_inference/outputs/LNCaP_scRNA_seacell_deg_results.pkl", "wb") as f:
    pickle.dump(all_degs, f)


# with open(f"{datadir}/TF_regulatory_inference/outputs/LNCaP_scRNA_seacell_deg_results.pkl", "rb") as f:
#     all_degs = pickle.load(f)

for gene_set in gene_sets:
    summary_df, all_gsea_results = seacells_compute_gsea(all_degs, gene_sets=gene_set)
    # Save summary
    summary_file = f"{datadir}/gsea_summary_{gene_set}_vs_dmso.csv"
    summary_df.to_csv(summary_file, index=False)
    # Save per-dataset GSEA results
    for dataset, gsea_result in all_gsea_results.items():
        gsea_result.res2d.to_csv(
            f"{datadir}/TF_regulatory_inference/outputs/gsea_results_{gene_set}_{dataset}_vs_dmso.csv"
        )

summary_df_hm = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/gsea_summary_MSigDB_Hallmark_2020_vs_dmso.csv")
terms = summary_df_hm['Pathway'].tolist()


auc_full = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/auc_scores_lncap_scRNA_msigdb_hallmark_20240307.csv")
auc_full_cols = auc_full.columns
common_rows = auc_full.index.intersection(rna.obs_names)

auc_full = auc_full[auc_full.index.isin(common_rows)]
auc_full.columns

def modify_term(term):
    # Replace special characters with '_', convert to uppercase, then remove duplicate underscores
    term = re.sub(r'\W', '_', term).upper()  # Replace non-word characters with '_'
    term = re.sub(r'_+', '_', term)          # Replace multiple underscores with a single '_'
    return term.strip('_')                   # Remove leading/trailing underscores if any

#  Apply the function to each term in the list
modified_terms = [modify_term(term) for term in terms]
modified_terms

filtered_columns = [col for col in auc_full.columns if col.replace("HALLMARK_", "") in modified_terms]

# Filter the 'auc_full' DataFrame to keep only columns of interest
filtered_auc_full = auc_full[filtered_columns]

auc_full_with_dataset = filtered_auc_full.copy()
auc_full_with_dataset['dataset'] = rna.obs['dataset'].values

# Group by dataset
mean_auc_per_dataset = auc_full_with_dataset.groupby('dataset').mean()

from matplotlib.colors import LinearSegmentedColormap
colors = ["#1a53ff", "#ffff99", "#cc0000"]
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

new_ticks = np.arange(0, 1.01, 0.25)  # Adjusting to range from 0 to 1 in steps of 0.25

# ------------------------------------------------------
# MSigDB Hallmark Heatmap (Supplementary Fig. 1G) 
# ------------------------------------------------------

# Create the clustermap 
g = sns.clustermap(mean_auc_per_dataset.T,
                   annot=False,
                   cmap=custom_cmap,  # Symmetric colormap (coolwarm)
                   fmt=".2f",
                   figsize=(8, 8),
                   cbar_kws={"label": "Score"},
                   standard_scale=0,
                   metric="euclidean",
                   method="average",
                   row_cluster=True,
                   col_cluster=True)


# Set new ticks on the color bar
g.cax.set_yticks(new_ticks)
g.cax.set_yticklabels([f"{t:.2f}" for t in new_ticks])
g.savefig(f"{datadir}/TF_regulatory_inference/outputs/dataset_msigdb_hallmark_heatmap_vs_dmso.pdf", dpi=300, bbox_inches="tight")
