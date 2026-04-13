import anndata as ad
import scanpy as sc
import os
import pandas as pd
import scglue
import numpy as np
from itertools import chain
from pyjaspar import jaspardb
import itertools
import pybedtools
import requests
import subprocess
import random
from pathlib import Path
datadir = Path.cwd() 

# download from http://download.gao-lab.org/GLUE/cisreg/JASPAR2022-hg38.bed.gz
motifbed = f"{datadir}/JASPAR2022-hg38.bed.gz"
date="20240307"
random.seed(321)
rna = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-"+date+".h5ad")
rna.var_names.name = "genes"
rna.var["name"] = rna.var_names
rna.obs_names.name = "cells"

# Fetch jaspar mofits
motif_bed = scglue.genomics.read_bed(motifbed)

# TFs matching the motifs present in the RNA data
tfs = pd.Index(motif_bed["name"]).intersection(rna.var_names)

# keeping the HVGs for target inference
genes = rna.var.query("highly_variable").index

# PySCENIC uses loom files as input
# create a loom with the HVGs and TFs and save the TFs to separate txt file
rna[:, np.union1d(genes, tfs)].write_loom(f"{datadir}/rna_{date}.loom")
np.savetxt(f"{datadir}/tfs_{date}.txt", tfs, fmt="%s")

print("building a coexpression-based draft network")

# command for building the coexpression networks:
# expects the arboreto_with_multiprocessing.py to be in the current dir 
# downloaded from https://github.com/aertslab/pySCENIC/

command = [
    "python","arboreto_with_multiprocessing.py", datadir+"/TF_regulatory_inference/outputs/rna_"+date+".loom",datadir+"/TF_regulatory_inference/outputs/tfs_"+date+".txt",
    "-o", datadir+"/TF_regulatory_inference/outputs/draft_grn_"+date+".csv", "-m", "grnboost2", "--seed", "0",
    "--num_workers", "20", "--cell_id_attribute", "cells",
    "--gene_attribute", "genes"
]

subprocess.run(command, check=True)

print("finished")
