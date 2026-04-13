
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import utils as ut
import json
import numpy as np
import random
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr
import matplotlib.cm as cm
import matplotlib.colors as colors
import os
import re
import warnings
warnings.filterwarnings("ignore")

random.seed(5325)
from pathlib import Path
datadir = Path.cwd()


def plot_correlations(df, y, gs_x,xlim=(-1.1, 1),ylim=(-1, 1),add_colorbar=True):
    # color scale for -1-1
    norm = colors.Normalize(vmin=-1, vmax=1)
    cmap = cm.get_cmap('seismic')
    n_cols = len(gs_x)
    fig, axes = plt.subplots(1, n_cols,figsize=(n_cols * 2, 2),squeeze=False)
    axes = axes[0]  # 1D array of axes
    df[y] = pd.to_numeric(df[y], errors='coerce')
    for j, x in enumerate(gs_x):
        ax = axes[j]
        df[x] = pd.to_numeric(df[x], errors='coerce')
        # Regression line
        sns.regplot(x=df[x], y=df[y],scatter=False, ax=ax,color='black',line_kws={'linewidth':2.5})
        # Spearman correlation
        corr, p = spearmanr(df[x], df[y])
        if p < 0.05:
            corr_color = cmap(norm(corr))
        else: 
            corr_color = "#b3b3b3" # grey if non-significant correlation
        # Scatterplot
        sns.scatterplot(x=df[x], y=df[y],s=30, alpha=0.8,color= corr_color,linewidth=0.2,legend=False,ax=ax, zorder =0)
        ax.set_title(f"R={corr:.2f}",fontsize=8)
        ax.set_xlabel(x, fontsize=6)
        ax.set_ylabel(y, fontsize=6)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    plt.tight_layout()
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, shrink=0.6, pad=0.02)
    cbar.set_label("Spearman R", fontsize=8)
    return fig



with open(f"{datadir}/data/processed_gene_sets.json", "r") as fp:
    gene_sets = json.load(fp)


cell_types = ['Basal', 'Duct luminal', 'Luminal', 'NE_identity']
ar_nepc_gs = ['AR_activity', 'Beltran_NEPC']
gs_y = ['SOX4_associated_resistance','SOX4_regulon','ARS_stress_linked_sensitive']

datasets = ['su2c', 'labrecque']

# -----------------------------------------
# Correlation scatter plots (Fig 3C and 3D) 
# -----------------------------------------

for dataset in datasets:
    df = pd.read_csv(f"{data_path}/bulkRNA_tumor/outputs/{dataset}_gsva_metadata.csv")
    df['Beltran_NEPC'] = df['Beltran_NEPC_UP'] - abs(df['Beltran_NEPC_DOWN'])
    filename = f"{data_path}/figures/{dataset}_corr_plots_SOX4_ARS.pdf"
    with PdfPages(filename) as pdf:
        for y in gs_y:
            #  Cell types page 
            fig1 = plot_correlations(df, y, cell_types)
            pdf.savefig(fig1, dpi=300)
            plt.close(fig1)
