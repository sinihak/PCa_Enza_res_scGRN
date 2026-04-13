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


main_dir = '/scratch/project_2014660/scGRN/scprinter'
work_dir = f'{main_dir}/seq2print'
frag_dir = f'{main_dir}/fragments'

random.seed(332218)

chromvar = anndata.read_h5ad(f'{work_dir}/chromvar_cisbp.h5ad')

# Perform dimension reduction of the single-cell motif scores object
sc.tl.pca(chromvar)
from cuml import UMAP
vec = UMAP(metric='cosine').fit_transform(chromvar.obsm["X_pca"])
chromvar.obsm['X_umap'] = vec

chromvar.write(f'{work_dir}/chromvar_cisbp.h5ad')

genome = scp.genome.hg38
printer = scp.load_printer(f'{work_dir}/lncap_scATAC_scprinter.h5ad', genome)

# Now fetch models:
import wandb
# Set your entity and project
entity = 'sini-hakkola-tampere-university'  # Replace with your W&B entity (username or team name)
project = 'scPrinter_seq_LNCaP_scATAC'  # Replace with your W&B project name

# Initialize the API
api = wandb.Api()

# Get the project
runs = api.runs(f"{entity}/{project}")

model_path = []
pretrain_models = []
for run in runs:
    if run.state != 'finished':
        continue

    if 'LNCaP' in run.tags and "Bulk" in run.tags:
        model_name = run.config['savename'] + '-' + run.name + '.pt'
        model_path.append(os.path.join(work_dir, "model", model_name))
        pretrain_models.append(model_name)
model_path

indir = '/scratch/project_2014660/scGRN/scprinter_inputs/'

seacells = pd.read_csv(os.path.join(indir, 'lncap_scATAC_seacells_formatted_for_scPrinter.csv'), index_col=0)
# barcode_groups = seacells.rename({'SEACell_alias':'group'}, axis=1)

cell_grouping, group_names = scp.utils.df2cell_grouping(printer, seacells)

cell_barcodes = np.array(chromvar.obs.index)

# run LoRA fine-tuning
lora_configs = []
embeddings = pd.DataFrame(chromvar.obsm["X_pca"], index=cell_barcodes)

for fold, model in enumerate(pretrain_models):
    lora_config = scp.tl.seq_lora_model_config(printer,
                                            region_path=f'{work_dir}/seq2print_merged_cleaned_narrowPeak.bed',
                                            cell_grouping=cell_grouping,
                                            group_names=group_names,
                                            embeddings=embeddings,
                                            genome=genome,
                                            pretrain_model=f'{work_dir}/model/{model}',
                                            overwrite_barcode=False,
                                            model_name=f'LNCaP_LoRA',
                                            fold=fold,
                                            model_config=f'{work_dir}/configs/LNCaP_fold{fold}.JSON',
                                            additional_lora_config={
                                            "lr":3e-5, # put smaller lr for further finetune, default: 3e-3
                                            "notes": "v3",
                                            "tags": ["LNCaP", "LoRA", f"fold{fold}"]},
                                            path_swap=(work_dir, ''),
                                            config_save_path=f'{work_dir}/configs/LNCaP_LoRA_fold{fold}.JSON')
    lora_configs.append(lora_config)

printer.close()

# This is to generate the terminal commands:
for fold in range(5):
    scp.tl.launch_seq2print(model_config_path=f'{work_dir}/configs/LNCaP_LoRA_fold{fold}.JSON',
                            temp_dir=f'{work_dir}/temp',
                            model_dir=f'{work_dir}/model',
                            data_dir=work_dir,
                            gpus=fold,
                            wandb_project='scPrinter_seq_LNCaP_scATAC',
                            verbose=False,
                            launch=False) # if launch=True, this command would launch the scripts directly,
                            # otherwise, it will just display the commands, which you should copy and run.