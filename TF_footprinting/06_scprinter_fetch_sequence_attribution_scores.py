import os, psutil
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

# here we fetch the sequence attribution scores for the de novo motif calling

main_dir = '/scratch/project_2014660/scGRN/scprinter'
work_dir = f'{main_dir}/seq2print'

random.seed(332218)
genome = scp.genome.hg38
printer = scp.load_printer(f'{work_dir}/lncap_scATAC_scprinter.h5ad', genome)

peaks = f'{work_dir}/regions.bed'

# Define the models to use
models = ['LNCaP_Bulk_fold0-crisp-pyramid-32.pt', 'LNCaP_Bulk_fold1-ancient-sea-33.pt','LNCaP_Bulk_fold2-likely-snow-34.pt','LNCaP_Bulk_fold3-smooth-paper-35.pt','LNCaP_Bulk_fold4-solar-mountain-36.pt']

model_path = [os.path.join(work_dir, "model", m) for m in models]

# generate count(accessibility)-based attribution scores

scp.tl.seq_attr_seq2print(
    genome=printer.genome,
    region_path=peaks,
    model_type='seq2print',
    model_path=model_path, # all folds
    gpus=[0,1,2,3],
    preset='count',
    overwrite=False,
    verbose=True,
    launch=False)


printer.close()