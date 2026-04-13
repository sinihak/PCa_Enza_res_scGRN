
import scprinter as scp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import time
import pandas as pd
import numpy as np
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


# here we print the commands that are later used for pretraining

## this will later be used as the base model for LoRA fine-tuning on pseudobulks/SEACells

# Specify directories we will use. Make sure you provide the full absolute path and not the relative path
main_dir = '/scratch/project_2014660/scGRN/scprinter'
work_dir = f'{main_dir}/seq2print'
frag_dir = f'{main_dir}/fragments'

random.seed(332218)

# Specify the reference genome
genome = scp.genome.hg38

printer = scp.load_printer(f'{work_dir}/lncap_scATAC_scprinter.h5ad', genome)

# Generate the config files for training seq2PRINT
import json
model_configs = []
if not os.path.exists(os.path.join(work_dir, 'configs')):
    os.makedirs(os.path.join(work_dir, 'configs'))
for fold in range(5):
    model_config= scp.tl.seq_model_config(printer,
                                     region_path=f'{work_dir}/seq2print_merged_cleaned_narrowPeak.bed',
                                     cell_grouping=printer.obs_names,
                                     group_names='LNCaP_Bulk',
                                     genome=printer.genome,
                                     fold=fold,
                                     overwrite_bigwig=False,
                                     model_name='LNCaP_Bulk',
                                     additional_config={
                                        "notes": "v3",
                                        "tags": ["LNCaP",
                                            "Bulk",
                                            f"fold{fold}"]},
                                     path_swap=(work_dir, ''),
                                     config_save_path=f'{work_dir}/configs/LNCaP_fold{fold}.JSON')
    model_configs.append(model_config)

printer.close()

for path in ['temp','model']:
    if not os.path.exists(os.path.join(work_dir, path)):
        os.makedirs(os.path.join(work_dir, path))

# Generate the commands for training the model
# 5-fold cross-validation
for fold in range(5):
    scp.tl.launch_seq2print(model_config_path=f'{work_dir}/configs/LNCaP_fold{fold}.JSON',
                            temp_dir=f'{work_dir}/temp',
                            model_dir=f'{work_dir}/model',
                            data_dir=work_dir,
                            gpus=fold, 
                            wandb_project='scPrinter_seq_LNCaP_scATAC', # wandb helps you manage loggins
                            verbose=False,
                            launch=False # if launch=True, this command would launch the scripts directly,
                            # otherwise, it will just display the commands, which you should copy and run.
                           )