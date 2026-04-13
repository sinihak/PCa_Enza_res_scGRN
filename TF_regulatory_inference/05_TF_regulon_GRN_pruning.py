import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
import pandas as pd
import numpy as np
from itertools import chain
from pyjaspar import jaspardb
import itertools
import seaborn as sns
import networkx as nx
from networkx.algorithms.bipartite import biadjacency_matrix
import pyscenic
import pyarrow
import pybedtools
import requests
from pyscenic.utils import load_motifs
import subprocess
import random
from pathlib import Path
datadir = Path.cwd() 


motifbed = f"{datadir}/JASPAR2022-hg38.bed.gz"

date="20240307"
random.seed(321)

# gene2peak guidance graph
gene2peak = nx.read_graphml(f"{datadir}/TF_regulatory_inference/outputs/guidance-0.1-qval-{date}.graphml.gz")

# the anndatas
rna = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-"+date+".h5ad")
atac = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/atac-emb-"+date+".h5ad")

genes = rna.var.query("highly_variable").index
peaks = atac.var.query("highly_variable").index

motif_bed = scglue.genomics.read_bed(motifbed)
tfs = pd.Index(motif_bed["name"]).intersection(rna.var_names)

#### TF cis-regulatory rankings bridged by the highly variable ATAC peaks ####

# ATAC-peak based ranking (ATAC regulatory element with TF motif)
# combines gene-peak and peak-TF connections into gene-TF cis-regulatory ranking

peak_bed = scglue.genomics.Bed(atac.var.loc[peaks])
peak2tf = scglue.genomics.window_graph(peak_bed, motif_bed, 0, right_sorted=True)
peak2tf = peak2tf.edge_subgraph(e for e in peak2tf.edges if e[1] in tfs)

# gene-tf ranking:
gene2tf_rank_glue = scglue.genomics.cis_regulatory_ranking(
    gene2peak, peak2tf, genes, peaks, tfs,
    # not needed as the peaks have fixed length
    # region_lens=atac.var.loc[peaks, "chromEnd"] - atac.var.loc[peaks, "chromStart"],
    random_state=0
)

#### TF cis-regulatory ranking with proximal promoters ####
# promoter-based ranking (if no regulatory ATAC peak)
# expanding to connect TSS-flanking regions (+/- 500 bp) with motif hits (genomic overlap)
# will be named based on the corresponding gene in the graph
flank_bed = scglue.genomics.Bed(rna.var.loc[genes]).strand_specific_start_site().expand(500, 500)
flank2tf = scglue.genomics.window_graph(flank_bed, motif_bed, 0, right_sorted=True)

# the gene-peak connections are replaced with gene-flank connections
# (self-loops, because each TSS flanking region has the same name as the corresponding gene)
gene2flank = nx.Graph([(g, g) for g in genes])
gene2tf_rank_supp = scglue.genomics.cis_regulatory_ranking(
    gene2flank, flank2tf, genes, genes, tfs,
    n_samples=0
)

# check five first: 
# gene2tf_rank_supp.iloc[:5, :5]

# _glue: ATAC-based, _supp: promoter based
gene2tf_rank_glue.columns = gene2tf_rank_glue.columns + "_glue"
gene2tf_rank_supp.columns = gene2tf_rank_supp.columns + "_supp"

# write the cis-regulatory rankings as feather files for pyScenic
scglue.genomics.write_scenic_feather(gene2tf_rank_glue, f"{datadir}/TF_regulatory_inference/outputs/glue_{date}.genes_vs_tracks.rankings.feather")
scglue.genomics.write_scenic_feather(gene2tf_rank_supp, f"{datadir}/TF_regulatory_inference/outputs/supp_{date}.genes_vs_tracks.rankings.feather")

# generate the annotation file
pd.concat([
    pd.DataFrame({
        "#motif_id": tfs + "_glue",
        "gene_name": tfs
    }),
    pd.DataFrame({
        "#motif_id": tfs + "_supp",
        "gene_name": tfs
    })
]).assign(
    motif_similarity_qvalue=0.0,
    orthologous_identity=1.0,
    description="placeholder"
).to_csv(f"{datadir}/TF_regulatory_inference/ctx_annotation_{date}.tsv"), sep="\t", index=False)

command = [
    "pyscenic", "ctx",
    datadir+"/TF_regulatory_inference/outputs/draft_grn_"+date+".csv",
    datadir+"/TF_regulatory_inference/outputs/glue_"+date+".genes_vs_tracks.rankings.feather",
    datadir+"/TF_regulatory_inference/outputs/supp_"+date+".genes_vs_tracks.rankings.feather",
    "--annotations_fname", datadir+"/TF_regulatory_inference/outputs/ctx_annotation_"+date+".tsv",
    "--expression_mtx_fname", datadir+"/TF_regulatory_inference/outputs/rna_"+date+".loom",
    "--output", datadir+"/outputs/TF_regulatory_inference/pruned_grn_"+date+".csv",
    "--rank_threshold", "400", # set based on the hvgs
    "--num_workers", "20",
    "--cell_id_attribute", "cells",
    "--gene_attribute", "genes"
]
subprocess.run(command, check=True)

# Load the pruned GRN data
tf2gene = load_motifs(f"{datadir}/TF_regulatory_inference/outputs/pruned_grn_" + date + ".csv")
# Get the target genes as a Series, keeping the transcription factor (TF) information
target_genes = tf2gene[("Enrichment", "TargetGenes")].apply(pd.Series).stack().reset_index(level="TF")

# Rename columns for clarity
target_genes.columns = ["TF", "TargetGene"]

# Keep only the unique target genes for each TF by grouping by "TF" and selecting unique target genes
unique_target_genes = target_genes.groupby("TF")["TargetGene"].unique().explode()

# Save the unique target genes for each TF to a CSV file
unique_target_genes.to_csv(f"{datadir}/TF_regulatory_inference/outputs/regulons_pruned_grn_" + date + ".csv", header=True)

# Get all unique transcription factors (TFs)
all_tfs = tf2gene.index.get_level_values("TF").unique()

# Print the number of transcription factors in the pruned network
print("Number of TFs in the pruned network:", len(all_tfs))
