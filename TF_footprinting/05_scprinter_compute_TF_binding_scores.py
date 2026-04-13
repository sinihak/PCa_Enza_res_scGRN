import sys
import os
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
import pyranges as pr
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
from scanpy.plotting.palettes import zeileis_28
from matplotlib.backends.backend_pdf import PdfPages
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
import warnings
warnings.filterwarnings("ignore")

random.seed(332218)
# Specify the reference genome. This must match that of your ATAC fragments file
genome = scp.genome.hg38

# Specify directories we will use. Make sure you provide the full absolute path and not the relative path
main_dir = '/scratch/project_2014660/scGRN/scprinter'
work_dir = f'{main_dir}/seq2print'

# Now fetch models:
import wandb
# Set your entity and project
entity = 'sini-hakkola-tampere-university'  # Replace with your W&B entity (username or team name)
project = 'scPrinter_seq_LNCaP_scATAC'  # Replace with your W&B project name

# Initialize the API
api = wandb.Api()

# Get the project
runs = api.runs(f"{entity}/{project}")


# Specify the reference genome. This must match that of your ATAC fragments file
genome = scp.genome.hg38

indir = '/scratch/project_2014660/scGRN/scprinter_inputs/'

# fetch the genes expressed in the scRNA-seq data (we are only interested in the TF binding scores of TFs that are expressed in the RNA data):
rna_genes = pd.read_csv(f'{indir}/lncap_scRNA_genes_20240307.csv')
rna_genes = list(rna_genes['x'])

# regions of interest were fetched like this:
# import scglue
# date = 
# rna = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-w-leiden-phases-tf-regulon-auc-seacells-aliases-{date}.h5ad")
## get genes, tss, and promoters
# TFs_to_keep = ["SOX4"]
# genes = scglue.genomics.Bed(rna.var.assign(name=rna.var_names))
# tss = genes.strand_specific_start_site()
# promoters = tss.expand(2000, 0)
# promoters = promoters[promoters.index.isin(TFs_to_keep)]
# roi = promoters.expand(150000,150000)
# roi.to_csv(f'{datadir}/outputs/TF_SOX4_regulatory_region_of_interest_150kb.bed',
#                   sep='\t', header=False, index=False)


# read the regions of interest

roi =  pr.read_bed(f'{indir}/TF_regulatory_inference/TF_SOX4_regulatory_region_of_interest_150kb.bed')

# fetch the peaks
merged_chromvar =  pr.read_bed(f'{work_dir}/chromvar_merged_cleaned_narrowPeak.bed')


overlapping = merged_chromvar.overlap(roi)
# convert to df
overlapping_df = overlapping.df

overlapping_df.to_csv(f'{work_dir}/scprinter_chromvar_peaks_overlapping_SOX4_regulatory_region_of_interest_150kb.csv',
                   sep='\t', header=False, index=False)

# read SEAcell/pseudobulk information
# The format should be index, barcode, group (i.e. seacell ID)
seacells = pd.read_csv(os.path.join(indir, 'lncap_scATAC_seacells_formatted_for_scPrinter.csv'), index_col=0)
selected_seacells = list(set(seacells['group']))

printer = scp.load_printer(f'{work_dir}/lncap_scATAC_scprinter.h5ad', genome)

cols = ['Chromosome', 'Start', 'End']
regions_df = overlapping_df[cols]

# save the regions overlapping the peaks
regions_df.to_csv(f'{work_dir}/SOX4_regulatory_regions_of_interest_150kb_peaks.bed',
                  sep='\t', header=False, index=False)

# convert to dictionary
regions_dict = {
    f"{row['Chromosome']}-{row['Start']}-{row['End']}": None
    for _, row in regions_df.iterrows()
}

# Define the models to use
lora_models = ['LNCaP_LoRA_fold0-fancy-dew-45.pt',
               'LNCaP_LoRA_fold1-hardy-eon-44.pt',
               'LNCaP_LoRA_fold2-drawn-durian-45.pt',
               'LNCaP_LoRA_fold3-bright-morning-45.pt',
               'LNCaP_LoRA_fold4-rare-butterfly-45.pt'
              ]

lora_model_path = [os.path.join(work_dir, "model", m) for m in lora_models]
import json

adata_tfbs = scp.tl.seq_tfbs_seq2print(seq_attr_count=None,
                      seq_attr_footprint=None, 
                      genome=printer.genome,
                      region_path=f'{work_dir}/SOX4_regulatory_region_of_interest_150kb_peaks.bed',
                      gpus=[0,1,2,3],
                      model_type='lora',
                      model_path=lora_model_path, # All folds/models
                      lora_config=json.load(open(f'{work_dir}/configs/LNCaP_LoRA_fold0.JSON', 'r')),
                      group_names=selected_seacells,
                      verbose=False,
                      launch=True,
                      return_adata=False, # To extract bigwig files for the next step
                      overwrite_seqattr=True,
                      save_key='LNCaP_LoRA_test', # and input a save_key
                      save_path=work_dir)
printer.close()

# Scan TF motifs across all regions to find motif matched sites

# Initialize motif set object
# bg = background nucleotide frequency, calculated from the given genome sequence
motifs = scp.motifs.JASPAR2022_core_Motifs(genome=genome, bg=[0.25] * 4)

# Prepare motif scanner. You can specify which TF motifs you want to scan using tf_genes. If tf_genes=None then use all motifs
motifs.prep_scanner()

regions = pd.read_csv(f'{work_dir}/SOX4_regulatory_region_of_interest_150kb_peaks.bed', sep='\t', header=None)

# Scan motif sites. This will return the exact genomic coordinates of motif matches
motif_sites = motifs.scan_motif(regions, verbose=True, clean=True)

# Reformat motif matches to a pandas DataFrame
motif_sites = pd.DataFrame(motif_sites)
motif_sites.iloc[:, 2] = motif_sites.iloc[:, 1] + motif_sites.iloc[:, 8]
motif_sites.iloc[:, 1] = motif_sites.iloc[:, 1] + motif_sites.iloc[:, 7]
motif_sites = motif_sites.iloc[:, [0,1,2,4]]
motif_sites.columns=["chrom", "start", "end", "TF"]
motif_sites.to_csv(f'{work_dir}/motif_sites_regulatory_regions_of_interest_150kb_peaks.csv')
motif_sitesquit

def fetch_bw(args):
    import pyBigWig as pw
    import numpy as np
    from tqdm import tqdm
    TFBS, bw, chrom_sizes = args
    chroms = np.array(TFBS['chrom'])
    starts = np.array(TFBS['start'])
    ends = np.array(TFBS['end'])
    res_all = {}
    with pw.open(bw, 'r') as f:
        bw_chroms = f.chroms()  # dict of chromosomes and their lengths in BigWig
        for chrom in tqdm(chrom_sizes):
            if chrom == 'chrY':
                continue
            if chrom not in bw_chroms: # skip chromosome if not present
                continue
            length = min(chrom_sizes[chrom], bw_chroms[chrom])
            res_all[chrom] = f.values(chrom, 0, length, numpy=True)
    vs = []
    for chr, left, right in zip(tqdm(chroms, mininterval=1), starts, ends):
        vs.append(np.nanmean(res_all[chr][left:right]))
    return vs

samples = selected_seacells
bigwig_dict = {sample:f"{work_dir}/{sample}_TFBS.bigwig" for sample in samples}
bigwig_dict
chrom_sizes = genome.chrom_sizes
args = [[motif_sites, bigwig_dict[sample], chrom_sizes] for sample in samples]

n_jobs = 9
import multiprocessing as mp
with mp.Pool(n_jobs) as pool:
    TFBS_scores = list(pool.imap(fetch_bw, args))
TFBS_scores = np.array(TFBS_scores).T
TFBS_scores = pd.DataFrame(TFBS_scores, columns=[f"TFBS_{sample}" for sample in samples])
TFBS_scores = pd.concat([motif_sites, TFBS_scores], axis=1)

# filter to keep only TF genes expressed in the RNA data
TFBS_scores = TFBS_scores[TFBS_scores['TF'].isin(rna_genes)]
TFBS_scores.to_csv(f'{work_dir}/plots/TFBS_scores_all_seacells_SOX4_regulatory_region_of_interest_150kb_all_folds.csv')