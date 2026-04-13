import scprinter as scp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import time
import pandas as pd
import numpy as np
import os
import random
import pickle
import torch
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
from scanpy.plotting.palettes import zeileis_28
from tqdm.contrib.concurrent import *
from tqdm.auto import *
import anndata
import scanpy as sc
import statistics as stat
import json
import csv
import re
import copy
import wandb
from sklearn.preprocessing import OneHotEncoder


# Specify directories we will use. Make sure you provide the full absolute path and not the relative path
main_dir = '/scratch/project_2014660/scGRN/scprinter'
work_dir = f'{main_dir}/seq2print'
frag_dir = f'{main_dir}/fragments'

random.seed(332218)

# let's construct the peak-by-cell count matrix using the peaks computed with the Chromvar preset
# First construct a peak-by-cell matrix of ATAC counts

genome = scp.genome.hg38
printer = scp.load_printer(f'{work_dir}/lncap_scATAC_scprinter.h5ad', genome)

# Fetch the cleaned peaks, save, it will be used in the next step
merged_chromvar = pd.DataFrame(printer.uns["peak_calling"]['merged_chromvar'][:])
merged_chromvar.to_csv(f'{work_dir}/regions.bed',
                     sep='\t', header=False, index=False)


adata = scp.pp.make_peak_matrix(printer,
                       regions=merged_chromvar,
                       region_width=300,
                       cell_grouping=None,
                       group_names=None,
                       sparse=True)

adata.write(f'{work_dir}/lncap_scATAC_cell_peak.h5ad')

printer.close()

# this anndata will be used for remocing the low coverage-regions for computing the TF binding scores later

# adata = anndata.read_h5ad(f'{work_dir}/cell_peak.h5ad')
coverage = adata.X.sum(axis=0)
adata = adata[:, coverage > 0]

# Calculate chromVAR motif scores using either GPU (device = "cuda", much faster) or CPU (device = "cpu", slower)

device = "cuda"

if device == "cuda":
    import warnings
    warnings.filterwarnings("ignore")
    import scanpy as sc
    import anndata
    import cupy as cp
    import cupyx as cpx
    import time
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    rmm.reinitialize(
        managed_memory=True, # Allows oversubscription
        pool_allocator=True, # default is False
        devices=0, # GPU device IDs to register. By default registers only GPU 0.
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)


# Sample background peaks for each peak
scp.chromvar.sample_bg_peaks(adata,
                             genome=genome,
                             method='chromvar',
                             niterations=250)

# Scan motifs
motif = scp.motifs.FigR_Human_Motifs(genome,
                                     bg=list(adata.uns['bg_freq']),
                                     n_jobs=100,
                                     pvalue=5e-5, mode='motifmatchr')
motif.prep_scanner(None, pvalue=5e-5)
motif.chromvar_scan(adata)

# Compute motif scores for single cells
chromvar = scp.chromvar.compute_deviations(adata, chunk_size=50000, device=device)

# Save for later use
chromvar.write(f'{work_dir}/chromvar_cisbp.h5ad')