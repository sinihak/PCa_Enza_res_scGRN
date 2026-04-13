
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import re

dataset_palette = {
    'dmso': '#1D61AE',
    'enz48h': '#E56D04',
    'resB': '#6C086D',
    'resA': '#367250'
}

phase_palette = {
    'G1': '#FF4D4D',
    'S': '#66B3FF',
    'G2M': '#00FF7F'
}


def load_regulons_from_csv(filename):
    import csv
    regulons = {}
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)

        # Skip the header
        next(reader)

        # Read each row and reconstruct the dictionary
        for row in reader:
            key = row[0]  # TF or group name is in the first column
            target_genes = row[1].split(',')  # Split the comma-separated string into a list
            regulons[key] = target_genes


def enrichment_analysis(datadir, all_genes, targets, padj_th,gene_sets, suffix=''):
    """
    Perform enrichment analysis for target genes in each TF group.

    Parameters:
    datadir (str): Directory where figures and results will be saved.
    all_genes (list): List of all possible genes.
    targets (dict): Dictionary where keys are TF group names and values are lists of target genes.
    Returns:
    None
    """
    # Ensure the output directories exist
    os.makedirs(os.path.join(datadir, 'figures'), exist_ok=True)
    os.makedirs(os.path.join(datadir, 'outputs'), exist_ok=True)
    # go_bp_2024 = '/scratch/spsiha/c5.go.bp.v2024.1.Hs.symbols.gmt'

    # Iterate over each TF group
    for group_name, gene_list in targets.items():
        # Ensure gene_list is a list of target genes
        if isinstance(gene_list, list):
            # Convert list to set to ensure unique target genes
            tf_group_targets = set(gene_list)

            # Perform enrichment analysis
            enr = gp.enrichr(
                gene_list=list(tf_group_targets),
                gene_sets=gene_sets,
                background=all_genes,
                outdir=None
            )

            # Filter significant results
            enr_sig = enr.results[enr.results['Adjusted P-value'] < padj_th]

            # Generate file paths
            if suffix:
                save_csv = os.path.join(datadir, f'outputs/gs_enrichment_{group_name}_{suffix}.tsv')
            else:
                save_csv = os.path.join(datadir, f'outputs/gs_enrichment_{group_name}.tsv')

            # save signicant terms to csv files and generate plots of top terms
            enr_sig.to_csv(save_csv, index=False, header=True, sep='\t')
        else:
            print(f"Warning: The value for '{group_name}' is not a list of genes. Skipping this TF group.")


def prepare_dotplot_data_with_top_terms(results, gene_sets, top_term, padj_th):
    """
    Prepare data for the dot plot based on the specified gene sets and filtering criteria.
    Parameters:
    results (pd.DataFrame): The concatenated enrichment results for all TF groups.
    gene_sets (list): List of gene sets to include in the dot plots.
    top_term (int): Number of top terms to include for each gene set.
    pval_th (float): Adjusted p-value threshold for filtering significant results.
    Returns:
    pd.DataFrame: A DataFrame prepared for plotting.
    """
    # Calculate -log10(Adjusted P-value)
    results["-log10(Adjusted P-value)"] = -np.log10(results["Adjusted P-value"])
    # Identify the top terms across all TF groups
    top_terms = get_top_terms(results, gene_sets, top_term)
    # Filter results to include only top terms across all groups
    filtered_results = results[(results["Term"].isin(top_terms)) & (results["Adjusted P-value"] < padj_th)].copy()
    # Calculate the number of contributing genes
    filtered_results.loc[:, "Gene Count"] = filtered_results["Genes"].apply(lambda x: len(x.split(";")))
    return filtered_results


def get_top_terms(df, gene_sets, top_term):
    """
    Identify the top terms across all TF groups for each gene set.
    Parameters:
    df (pd.DataFrame): The concatenated enrichment results for all TF groups.
    gene_sets (list): List of gene sets to include in the dot plots.
    top_term (int): Number of top terms to include for each gene set.
    Returns:
    set: A set of unique top terms across all TF groups.
    """
    top_terms = set()
    for gene_set in gene_sets:
        gene_set_df = df[df["Gene_set"] == gene_set]
        for group in df["TF Group"].unique():
            group_df = gene_set_df[gene_set_df["TF Group"] == group]
            # Use the nlargest function to get the top terms based on -log10(Adjusted P-value)
            top_terms.update(group_df.nlargest(top_term, "-log10(Adjusted P-value)")["Term"])
    return top_terms

def plot_dotplot(df, gene_set, datadir, top_term):
    """
    Create and save a dot plot for the enrichment results with dynamic dot size annotation.
    Parameters:
    df (pd.DataFrame): The prepared dot plot data.
    gene_set (str): The gene set to plot.
    datadir (str): Directory to save the plot.
    """
   # Filter data for the specific gene set
    plot_data = df[df["Gene_set"] == gene_set]
    # Dynamically calculate the plot height
    plot_height = max(3, 0.3 * len(plot_data["Term"].unique()))
    # Create the dot plot
    plt.figure(figsize=(10, plot_height))
    ax = sns.scatterplot(
        data=plot_data,
        x="TF Group",
        y="Term",
        size="Gene Count",  # Size of the dots based on the number of contributing genes
        hue="-log10(Adjusted P-value)",  # Color of the dots based on -log10(Adjusted P-value)
        palette="inferno",
        sizes=(20, 300),
        legend = True,
        edgecolor="black"
    )
    # Add a custom legend for dot sizes
    min_gene_count = plot_data["Gene Count"].min()
    max_gene_count = plot_data["Gene Count"].max()
    # Generate 5 evenly spaced values for size ticks
    size_ticks = np.linspace(min_gene_count, max_gene_count, 6)  # 5 evenly spaced gene counts
    size_mapping = [20 + (300 - 20) * ((val - min_gene_count) / (max_gene_count - min_gene_count)) for val in size_ticks]
    # Create scatter points for the size scale legend
    for size, val in zip(size_mapping, size_ticks):
        ax.scatter([], [], s=size, color="black", edgecolor="black", label=f"{int(val)} Genes")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)
    # Customize the plot
    ax.set_title(f"{gene_set}", fontsize=13)
    plt.xticks(rotation=45, ha="right")
    # Adjust space around the plot
    plt.subplots_adjust(left=0.1, right=0.35, top=0.9, bottom=0.1)
    # Save the plot
    save_path = f"{datadir}/figures/dotplot_{gene_set}_top_{top_term}_terms.pdf")
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()