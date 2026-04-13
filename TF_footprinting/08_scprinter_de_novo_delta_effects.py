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

main_dir = '/scratch/project_2014660/scGRN/scprinter'
work_dir = f'{main_dir}/seq2print'

random.seed(332218)
genome = scp.genome.hg38

peaks = f'{work_dir}/regions.bed'


models = ['LNCaP_Bulk_fold0-crisp-pyramid-32.pt',
'LNCaP_Bulk_fold1-ancient-sea-33.pt',
'LNCaP_Bulk_fold2-likely-snow-34.pt',
'LNCaP_Bulk_fold3-smooth-paper-35.pt',
'LNCaP_Bulk_fold4-solar-mountain-36.pt'
]

model_path = [os.path.join(work_dir, "model", m) for m in models]

# generate the delta effects for the de novo motifs

de = scp.tl.delta_effects_seq2print(model_path,
                               genome=scp.genome.hg38,
                               region_path=peaks,
                               motifs= f'{work_dir}/modisco/modisco.count.h5', # path to the de novo motifs
                               prefix='count',
			                   gpus=[0,1,2,3],
                               launch = False,
                               save_path=f'{work_dir}/modisco/delta_effects_count',
                               plot=True, vmin='auto', vmax='auto') # 'auto': 5th percentile is used

