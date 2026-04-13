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
print(scp.__version__)
random.seed(332218)

genome = scp.genome.hg38

main_dir = '/scratch/project_2014660/scGRN/scprinter'
work_dir = f'{main_dir}/seq2print'


# create the modisco report of the de novo motifs
scp.tl.modisco_report(modisco_h5=f'{work_dir}/modisco/modisco.count.h5',
                     save_path=f'{work_dir}/modisco_count',
                     meme_motif=scp.datasets.FigR_motifs_human_meme,
                     top_n_matches=3,
                     delta_effect_path=f'{work_dir}/modisco/delta_effects_count',
                     motif_prefix='count')

# define the models to use
models = ['LNCaP_Bulk_fold0-crisp-pyramid-32.pt',
'LNCaP_Bulk_fold1-ancient-sea-33.pt',
'LNCaP_Bulk_fold2-likely-snow-34.pt',
'LNCaP_Bulk_fold3-smooth-paper-35.pt',
'LNCaP_Bulk_fold4-solar-mountain-36.pt'
]

model_path = [os.path.join(work_dir, "model", m) for m in models]

peaks = f'{work_dir}/regions.bed'

print('Calling de novo hits for counts')
hits_count = scp.tl.seq_denovo_callhits(modisco_output=f'{work_dir}/modisco/modisco.count.h5',
                           model_path=model_path,
                           region_path=peaks,
                           device='cuda:0',
                           preset='count',
                           save_path=f'{work_dir}/modisco/finemo_count',
                           overwrite=False,
                           verbose=True,
                           launch=True,
                           return_hits=True)

hits_count.to_csv(f"{work_dir}/modisco/finemo_count/de_novo_hits_count.tsv",sep="\t", index=False)

# Set the GPU device for calculating the motif scores for the de novo motifs

import cupy as cp
import cupyx as cpx
import time
import warnings
warnings.filterwarnings("ignore")
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
rmm.reinitialize(
    managed_memory=True, # Allows oversubscription
    pool_allocator=True, # default is False
    devices=0, # GPU device IDs to register. By default registers only GPU 0.
)

cp.cuda.set_allocator(rmm_cupy_allocator)

# run ChromVAR for the de novo hits
layers= ['count']
all_hits = [hits_count]

adata = anndata.read_h5ad(f'{work_dir}/lncap_scATAC_cell_peak.h5ad')

for hits, layer in zip(all_hits, layers):
    print(f'starting to run Chromvar on {layer}')
    # work on the copy of the main object
    ad = adata.copy()
    # Take the hits, and create a motif matching data matrix and motif name list for chromVAR
    motif_uniq = np.sort(hits['motif_name'].unique())
    motif2id = {m:i for i,m in enumerate(motif_uniq)}
    ids = [motif2id[m] for m in hits['motif_name']]
    match_mm = np.zeros((ad.shape[1], len(motif_uniq)))
    match_mm[hits['peak_id'], ids] += 1
    ad.varm['motif_match'] = match_mm
    motif_uniq = [f'{layer}_{xx}' for xx in motif_uniq]
    ad.uns['motif_name'] = motif_uniq
    # Filter the low coverage regions (after adding the hits, because all peaks were used for computing them)
    coverage = ad.X.sum(axis=0)
    adata = ad[:, coverage > 0]
    # background peaks
    scp.chromvar.sample_bg_peaks(ad,
                             genome=scp.genome.hg38,
                             method='chromvar',
                             niterations=250)
    # now run chromVAR
    chromvar_denovo = scp.chromvar.compute_deviations(ad, chunk_size=50000, device='cuda')
    print('saving the chromvar scores')
    chromvar_denovo.write(f'{work_dir}/chromvar_de_novo_{layer}.h5ad')