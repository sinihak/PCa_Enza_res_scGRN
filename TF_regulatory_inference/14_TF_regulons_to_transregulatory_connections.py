
import os
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.stats import fisher_exact
import statsmodels.stats.contingency_tables as smc
import utils as ut
from statsmodels.stats.multitest import multipletests
import re
from pathlib import Path
datadir = Path.cwd() 

date= "20240307"
num="250"
GRNqval="0.1"

# parameters for DE testing
pval_cutoff=0.05
log2fc_min=0.5

### first, let"s calculate TF regulon DE genes ####

rna = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/lncap_scRNA_unscaled_{date}.h5ad")

# Load target genes
target_genes_df = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/regulons_pruned_grn_{date}.csv")
target_genes_df['genes'] = target_genes_df.iloc[:, 1].apply(lambda x: re.findall(r"'(.*?)'", x))
# Create regulons dictionary
regulons = target_genes_df.groupby("TF")["genes"].apply(list).to_dict()
TFs_to_keep = ["ZNF682","E2F8","E2F7","E2F1","E2F2", "YY1", "NFYB", "ZNF669", "ZNF574", "MSANTD3",
"CEBPG", "ATF4", "CREB3", "ZIC2", "IRF1", "SREBF1", "HIF1A", "SOX4", "FOXO3", "MYC","SREBF2", "ZNF282", "RREB1"]


# Filter regulons to keep only the specified TFs
regulons = {k: regulons[k] for k in TFs_to_keep if k in regulons}
# Extract unique genes
genes = list(set(gene for genes_list in regulons.values() for genes in genes_list for gene in genes))

flattened_regulons = {TF: [gene for sublist in genes for gene in sublist] for TF, genes in regulons.items()}


comparisons = [
    ("resA", "enz48h"),
    ("resA", "dmso"),
    ("resB", "enz48h"),
    ("resB", "dmso"),
    ("enz48h", "dmso"),
    ("dmso", "rest")
]

keys_to_check = ["stress_responses", "ARS_associated", "resistant", "cell_cycle_related"]

# Define the functions:

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
    sc.tl.rank_genes_groups(rna_temp, groupby=groupby, groups=[group], reference=reference, method="wilcoxon")
    df = sc.get.rank_genes_groups_df(rna_temp, group=group, pval_cutoff=pval_cutoff)
    df = df.dropna(subset=["logfoldchanges"])
    del rna_temp
    return df

# this is to calculate the DEGs for all comparisons defined above
def compute_all_DEGs(rna, comparisons, pval_cutoff):
    """Compute DEGs for all comparisons and return a dictionary."""
    deg_results = {}
    for group, ref in comparisons:
        key = f"{group}_vs_{ref}"
        df = compute_DEGs(rna, groupby="dataset", group=group, reference=ref, pval_cutoff=pval_cutoff)
        df["comparison"] = key
        deg_results[key] = df
    return deg_results

def filter_up_down(df, log2fc_th=log2fc_min):
    """Return upregulated and downregulated DEGs."""
    up = df[df["logfoldchanges"] > log2fc_th]
    down = df[df["logfoldchanges"] < -log2fc_th]
    return up, down


def filter_regulons(regulons, gene_set):
    """Filter regulons to include only genes present in gene_set."""
    return {reg: list(set(g for genes in gene_lists for g in genes if g in gene_set))
            for reg, gene_lists in regulons.items() if any(g in gene_set for genes in gene_lists for g in genes)}

def filter_regulon_groups(gene_set, regulon_groups, keys):
    """Filter regulon groups by set of genes (eg. DEGs)."""
    filtered = {}
    for key in keys:
        if key in regulon_groups:
            genes_lists = regulon_groups[key]
            filtered_genes = set(gene for gene in genes_lists if gene in gene_set)
            if filtered_genes:
                filtered[key] = list(filtered_genes)
    return filtered

################# 
# MAIN WORKFLOW #
#################

# Compute DEGs

deg_results = compute_all_DEGs(rna, comparisons, pval_cutoff)

# Save combined DEGs
dfs = pd.concat(deg_results.values(), axis=0)
dfs.to_csv(f"{datadir}/TF_regulatory_inference/lncap_scRNA_dataset_de_genes.csv")

# Filter upregulated genes for each comparison
# filter_up_down(df)[0] returns only the first output, i.e. the upregulated genes
upregulated = {key: filter_up_down(df)[0] for key, df in deg_results.items()}


# Extract the DE genes for the resistant conditions
resA_genes = {gene for comp in ["resA_vs_enz48h", "resA_vs_dmso"] for gene in upregulated[comp]["names"]}
resB_genes = {gene for comp in ["resB_vs_enz48h", "resB_vs_dmso"] for gene in upregulated[comp]["names"]}

# ------------------------------------------------------------
# TF regulon group association with RES-A/B DE genes (Fig. 1I) 
# ------------------------------------------------------------

groups = ["stress_responses", "cell_cycle_related","resistant", "ARS_associated"]
stats = []
for group in groups:
    group_genes = set(regulon_groups[group])
    resA_set = set(resA_genes) & group_genes # DE genes in ResA that are in regulon group
    resB_set = set(resB_genes) & group_genes # DE genes in ResB that are regulon group
    # all genes in group
    # Fisher"s exact test comparing DE CC proportions between ResA and ResB
    resA_non = len(set(resA_genes) - group_genes)    
    resB_non = len(set(resB_genes) - group_genes)
    contingency_table = [[len(resA_set), resA_non], [len(resB_set), resB_non]]
    # Calculate statistics
    odds_ratio, p_value = fisher_exact(contingency_table)
    table = smc.Table2x2(contingency_table)
    # compute confidence intervals
    ci_low, ci_high = table.oddsratio_confint(method="exact",alpha=0.05)
    stats.append({"group": group, "odds_ratio": odds_ratio, "p_value": p_value, "ci_low":ci_low, "ci_high":ci_high})

stats_df = pd.DataFrame(stats)

# Sort groups 
stats_df = stats_df.sort_values("p_value", ascending=False)
# Create figure
plt.figure(figsize=(4.5, 3))
# Y positions
y_pos = range(len(stats_df))
# Plot points and CIs
plt.errorbar(
    stats_df["odds_ratio"],
    y_pos,
    xerr=[
        stats_df["odds_ratio"] - stats_df["ci_low"],
        stats_df["ci_high"] - stats_df["odds_ratio"]
    ],
    fmt="o",
    capsize=4,
    markersize=6,
    ecolor="#2C1440",
    linestyle="none",
    markerfacecolor="#2C1440",
    markeredgecolor="#2C1440"
)
# Add y-axis labels
plt.yticks(y_pos, stats_df["group"])
# Add vertical reference line at OR = 1
plt.axvline(1, color="grey", linestyle="--")
plt.xlabel("Odds ratio")
plt.tight_layout()
plt.savefig(f"{datadir}/TF_regulatory_inference/figures/LNCaP_scRNA_regulon_groups_resA_resB_forest_plot.pdf"), dpi=300, pad_inches=0.5)
plt.close()


# load the TF regulons and coexpression modules
tf_regulons = ut.load_regulons_from_csv(f"{datadir}/TF_regulatory_inference/outputs/TF_regulons.csv")

MEs = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_scRNA_module_eigengenes.csv")
MEs.index = MEs["Unnamed: 0"]
MEs = MEs.drop("Unnamed: 0", axis=1)
# module_names = MEs.columns

rna.obs = rna.obs.join(MEs)

# ----------------------------------------------------------
# TF regulon association with coexpression modules (Fig. 1J) 
# ---------------------------------------------------------

modules = pd.read_csv(f"{datadir}/TF_regulatory_inference/outputs/coexpression_modules.csv")

# the entire gene universe will be used as they all had a potential to be included in the coexpression modules
all_genes = set(rna.var.index)

# Prepare results storage
results = []

# Iterate over modules and regulons
for module in modules["module"].unique():
    # skip the "grey" module
    if module == "grey":
        continue
    module_genes = set(modules[modules["module"] == module]["gene_name"])
    for tf, regulon_genes in tf_regulons.items():
        regulon_genes = set(regulon_genes)  # Convert to set
        # Compute overlap
        overlapping_genes = module_genes & regulon_genes
        a = len(overlapping_genes)  # Overlapping genes
        b = len(regulon_genes - module_genes)  # In regulon, not in module
        c = len(module_genes - regulon_genes)  # In module, not in regulon
        d = len(all_genes - (module_genes | regulon_genes))  # Neither in module nor regulon
        # Fisher’s Exact Test
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        # Store result
        results.append({
            "TF": tf,
            "Module": module,
            "Overlap": a,
            "P-value": p_value,
            "Overlapping Genes": ", ".join(overlapping_genes)  # Store overlapping genes as a comma-separated string
        })

# Convert results to a DataFrame
results_df = pd.DataFrame(results)
# Multiple testing correction (Benjamini-Hochberg)
results_df["Adjusted P-value"] = multipletests(results_df["P-value"], method="fdr_bh")[1]
# Save to CSV
results_df.to_csv(f"{datadir}/TF_regulatory_inference/outputs/TF_regulon_module_enrichment_results.csv"), index=False)
# Display significant results

sig_results = results_df[results_df["Adjusted P-value"] < 0.05]

results_df_sorted = sig_results.sort_values(by="Overlap", ascending=False)
pivot_df = results_df_sorted.pivot(index="TF", columns="Module", values="Overlap").fillna(0)

# Sort TFs based on the total overlap (sum across modules)
pivot_df["Total_Overlap"] = pivot_df.sum(axis=1)
pivot_df = pivot_df.sort_values(by="Total_Overlap", ascending=False)
pivot_df = pivot_df.drop(columns=["Total_Overlap"])  # Drop the helper column

# fetch pvals

# Plot stacked bars
ax = pivot_df.plot(kind="barh", stacked=True, figsize=(4, 3), cmap="tab20_r")
plt.xlabel("# of overlapping genes")
plt.ylabel("TF regulon")
plt.legend(title="Module", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.gca().invert_yaxis()  # Ensure highest overlap is at the top
plt.savefig(f"{datadir}/TF_regulatory_inference/figures/TF_regulon_overlap_with_gene_modules.pdf", dpi=300, bbox_inches="tight")


### SOX4 regulon overlaps with sc-M9 and MYC, SREBF2 and ATF4 with sc-M1
### these and regulon DE genes are used to generate gene sets to represent
### trans-regulatory connections, that will be tested with clinical datasets


# Strict threshold for SOX4 targets
# Most altered genes (log2FC >= 1)
deg_resistant_most_altered = {gene for name, df in upregulated.items() if "res" in name for gene in df.loc[df["logfoldchanges"] >= 1, "names"]}
most_altered_regulons = filter_regulons(regulons, deg_resistant_most_altered)
# SOX4 most altered regulon genes

print(most_altered_regulons["SOX4"])
# "MBNL2", "ZMAT1", "TUBA1A", "MET", "ZNF608", "KCNMB4", "CD24", "MARCKS", "ZSWIM5", "SESN3", "C16orf89", "TMSB4X" 


# We now have newly defined SOX4 regulon, only containing the DE genes: 
sox4_regulon = ["MBNL2", "ZMAT1", "TUBA1A", "MET", "ZNF608", "KCNMB4", "CD24", "MARCKS", "ZSWIM5", "SESN3", "C16orf89", "TMSB4X"]


m9 = modules[modules['Module']=='sc-M9']
m9_genes = list(m9['Genes'])
m9_genes = m9_genes[0].split(",")

# SOX4-linked resistance
SOX4_linked_resistance = list(set(m9_genes + sox4_regulon))



deg_dmso = {gene for gene in upregulated['dmso_vs_rest']["names"]}
deg_enz  = {gene for gene in upregulated['enz48h_vs_dmso']["names"]}

m1 = modules[modules['Module']=='sc-M1']
m1_genes = list(m1['Genes'])
m1_genes = m1_genes[0].split(",")

# Load regulon groups and filter
regulon_groups = ut.load_regulons_from_csv(f"{datadir}/TF_regulatory_inference/outputs/TF_regulons_groups.csv")

# filter regulon groups to contain only the "sensitive" DEGs
filtered_regulon_groups_dmso = filter_regulon_groups(deg_dmso, regulon_groups, keys_to_check)
filtered_regulon_groups_enz = filter_regulon_groups(deg_enz, regulon_groups, keys_to_check)

# Merge sensitive regulon groups
filtered_regulon_groups_sens = {}

for key in set(filtered_regulon_groups_dmso) | set(filtered_regulon_groups_enz):
    filtered_regulon_groups_sens[key] = list(
        set(filtered_regulon_groups_dmso.get(key, []) + filtered_regulon_groups_enz.get(key, []))
)


ARS_stress_linked_sensitive = list(set(m1_genes + filtered_regulon_groups_sens["ARS_associated"] + filtered_regulon_groups_sens["stress_responses"]))

# ----------------------------------------------------------------------------
# Visualize selected coexpression modules as UMAP (Supplementary Fig 1K and 1L) 
# ----------------------------------------------------------------------------

n_modules = len([m for m in module_names if m != "grey"])  # Exclude "grey"
n_cols = 2  # Number of columns in the grid
n_rows = int(np.ceil(n_modules / n_cols))  # Auto-calculate rows

# define some colors
colors = ["lightgray","darkgray", "gold", "orangered"]
from matplotlib.colors import LinearSegmentedColormap
# Create a linear segmented colormap
cmap = LinearSegmentedColormap.from_list("custom_palette", colors, N=50)

module_names = ["sc-M1", "sc-M9"] # these we will next plot as UMAPs

# Define grid layout (adjust rows/cols based on number of modules)
output_pdf = f"{datadir}/TF_regulatory_inference/figures/lncap_scRNA_module_eigengenes_subset.pdf"
with PdfPages(output_pdf) as pdf:
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 3))  # Adjust figure size
    axes = axes.flatten()  # Flatten for easy indexing
    plot_idx = 0
    for module in module_names:
        if module == "grey":
            continue  # Skip the "grey" module
        # Ensure the module data is numeric
        rna.obs[module] = pd.to_numeric(rna.obs[module], errors="coerce")
        # Check if there are valid numeric values
        if rna.obs[module].isna().all():
            print(f"Skipping {module}: No valid numeric data.")
            continue
        # Compute symmetric limits
        upper = rna.obs[module].mean() + 3 * rna.obs[module].std() # this approach is borrowed from hdWGCNA
        lower = rna.obs[module].mean() - 3 * rna.obs[module].std()
        abs_min = min(abs(lower), abs(upper))
        upper_limit = abs_min
        lower_limit = -abs_min
        # Plot each module on a separate subplot
        sc.pl.embedding(rna, color=module, basis="X_umap_comb", size=15,
                        vmin=round(lower_limit), vmax=round(upper_limit),
                        cmap=cmap, ax=axes[plot_idx], show=False)
        axes[plot_idx].set_title(module)  # Add title per subplot
        plot_idx += 1
    # Hide unused subplots if modules don't fill the grid
    for ax in axes[plot_idx:]:
        ax.axis("off")
    pdf.savefig(fig, bbox_inches="tight", dpi=300)
    plt.close(fig)  # Close to free memory
