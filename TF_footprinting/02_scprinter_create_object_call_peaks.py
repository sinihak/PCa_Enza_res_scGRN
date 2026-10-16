import os
os.environ["LOKY_USE_SEMAPHORE"] = "0"  # or "False"
import sys
import random
import scprinter as scp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import time
import pandas as pd
import numpy as np
import pickle
import random
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
from sklearn.preprocessing import OneHotEncoder
from joblib.externals.loky import get_reusable_executor
import tempfile
import warnings

warnings.filterwarnings("ignore")

random.seed(332218)

import os
import multiprocessing as mp

# Also set wandb directory if used
os.environ["WANDB_DIR"] = "/scratch/project_2014660/wandb"
os.makedirs(os.environ["WANDB_DIR"], exist_ok=True)

# Specify the reference genome (should match the fragments files)
genome = scp.genome.hg38

indir='/scratch/project_2014660/scGRN/scprinter_inputs/'

dmso_frags = os.path.join(indir, 'Sample1/fragments.tsv.gz')
enz48h_frags = os.path.join(indir, 'Sample2/fragments.tsv.gz')
resA_frags = os.path.join(indir, 'Sample3/fragments.tsv.gz')
resB_frags = os.path.join(indir, 'Sample4/fragments.tsv.gz')
fragment_files = [dmso_frags, enz48h_frags, resA_frags, resB_frags]

# Use the barcodes that have passed the QC filters:
atac = sc.read_h5ad(os.path.join(indir, 'atac-emb-w-leiden-phases-seacells-aliases-20240307.h5ad'))

import re

def strip_prefix(barcode_index, sample_prefix):
    return barcode_index[barcode_index.str.startswith(sample_prefix)].str.replace(f'^{sample_prefix}_', '', regex=True)

dmso_barcodes = strip_prefix(atac.obs.index, 'dmso')
enz48h_barcodes = strip_prefix(atac.obs.index, 'enz48h')
resA_barcodes = strip_prefix(atac.obs.index, 'resA')
resB_barcodes = strip_prefix(atac.obs.index, 'resB')


barcodes = [dmso_barcodes.tolist(), enz48h_barcodes.tolist(), resA_barcodes.tolist(), resB_barcodes.tolist()]

main_dir = '/scratch/project_2014660/scGRN/scprinter'
work_dir = f'{main_dir}/seq2print'
if not os.path.exists(work_dir):
    os.system("mkdir -p " + work_dir)
frag_dir = f'{main_dir}/fragments'
if not os.path.exists(frag_dir):
    os.system("mkdir -p " + frag_dir)

samples = ['dmso', 'enz48h', 'resA', 'resB']

# Initialize the scPrinter object
# When you finish using the object, run printer.close() otherwise you won't be able to load it properly next time.

import time
start = time.time()
printer_path = os.path.join(work_dir, 'lncap_scATAC_scprinter.h5ad')


def main():
    if os.path.exists(printer_path):
        printer = scp.load_printer(f'{work_dir}/lncap_scATAC_scprinter.h5ad', genome)
    else:
        printer = scp.pp.import_fragments(
            path_to_frags=fragment_files,
            barcodes=barcodes,
            sample_names=samples,
            savename=printer_path,
            genome=genome,
            sorted_by_barcode=False,
            low_memory=True,
            n_jobs=1
        )
    scp.pp.call_peaks(
        printer=printer,
        frag_file=fragment_files,
        cell_grouping=barcodes,
        group_names=[f"{s}_seq2print" for s in samples],
        iterative_peak_merging=True,
        merged_key="merged_seq2print",  # custom key
        merge_across_groups=True,
        preset='seq2PRINT',
        n_jobs=1
    )

    merged_seq2print = pd.DataFrame(printer.uns["peak_calling"]['merged_seq2print'][:])
    merged_seq2print.to_csv(f'{work_dir}/seq2print_merged_cleaned_narrowPeak.bed',
                 sep='\t', header=False, index=False)

    scp.pp.call_peaks(
        printer=printer,
        frag_file=fragment_files,
        cell_grouping=barcodes,
        group_names=[f"{s}_chromvar" for s in samples],
        iterative_peak_merging=True,
        merge_across_groups=True,
        preset='chromvar',
        n_jobs=1,
        merged_key="merged_chromvar",  # custom key
        overwrite=False
    )

    merged_chromvar = pd.DataFrame(printer.uns["peak_calling"]['merged_chromvar'][:])
    merged_chromvar.to_csv(f'{work_dir}/chromvar_merged_cleaned_narrowPeak.bed',
                 sep='\t', header=False, index=False)

    printer.close()

if __name__ == "__main__":
     main()