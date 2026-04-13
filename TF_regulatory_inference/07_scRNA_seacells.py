import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import SEACells
import matplotlib
import scipy
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.io import mmread
from matplotlib.backends.backend_pdf import PdfPages
import random
from pathlib import Path

# fetches current wd
datadir = Path.cwd() 

date="20240307"

rna = ad.read_h5ad(f"{datadir}/TF_regulatory_inference/anndata/rna-emb-w-leiden-phases-{date}.h5ad")

num = 250 ## approx. size of metacells (nCells/metacell)
nMetaCells = round(len(rna.obs.index)/num)
print("Number of metacells:", nMetaCells)
## Core parameters 
n_SEACells = nMetaCells ## ~X single cells for each metacell
build_kernel_on = "X_pca"

## Additional parameters
n_waypoint_eigs = 7  # number of eigenvalues to consider when initializing metacells

## Make model  
print("making the model")
model = SEACells.core.SEACells(rna, 
                      build_kernel_on=build_kernel_on, 
                      n_SEACells=n_SEACells,
                      n_neighbors=15,   
                      n_waypoint_eigs=n_waypoint_eigs,
                      convergence_epsilon = 1e-5)
model.construct_kernel_matrix()
M = model.kernel_matrix

## Initialize archetypes
model.initialize_archetypes()
print("fitting the model")
model.fit(min_iter=10, max_iter=70)


# Some diagnostics plots
SEACell_purity = SEACells.evaluate.compute_celltype_purity(rna, "leiden")
with PdfPages(f"{datadir}/TF_regulatory_inference/figures/leiden_purity_num_{num}_{date}.pdf") as pdf:
    # Plot the boxplot and save it to the PDF
    plt.figure(figsize=(4, 4))
    sns.boxplot(data=SEACell_purity, y="leiden_purity")
    plt.title("Cluster Purity")
    sns.despine()
    pdf.savefig(bbox_inches="tight")
    plt.close()


# Lower values of compactness suggest more compact/lower variance metacells.
compactness = SEACells.evaluate.compactness(rna, build_kernel_on)
compactness.to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_rna_compactness_num_{num}_{date}.csv")

with PdfPages(f"{datadir}/TF_regulatory_inference/metacell_compactness_num_{num}_{date}.pdf") as pdf:
    plt.figure(figsize=(4, 4))
    sns.boxplot(data=compactness, y="compactness")
    plt.title("Compactness")
    sns.despine()
    pdf.savefig(bbox_inches="tight")
    plt.close()

# Higher values of separation suggest better distinction between metacells.
separation = SEACells.evaluate.separation(rna, "X_pca", nth_nbr=1)
separation.to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_rna_separation_num_{num}_{date}.csv")
with PdfPages(f"{datadir}/TF_regulatory_inference/figures/metacell_separation_num_{num}_{date}.pdf") as pdf:
    plt.figure(figsize=(4, 4))
    sns.boxplot(data=separation, y="separation")
    plt.title("Separation")
    sns.despine()
    pdf.savefig(bbox_inches="tight")
    plt.close()
print("starting to plot diagnostics")


model_path = f"{datadir}//F_regulatory_inference/outputs/lncap_rna_seacells_model_num_{num}_{date}.pickle"
with open(model_path, "wb") as file:
    pickle.dump(model, file)


# replace the SEACell label with the archetype cell name 
soft_labels, weights = model.get_soft_assignments()
# The archetype label is in the first column
seacell_labels = soft_labels.iloc[:, 0]
# This assumes that the indices match, you can check that to be sure.

indices_match = soft_labels.index.equals(rna.obs.index)

if indices_match:
    print("The indices match.")
else:
    print("The indices do not match.")


rna.obs["SEACell"] = seacell_labels 

print("making the feature-metacell matrix")
SEACell_rna = SEACells.core.summarize_by_SEACell(rna, SEACells_label="SEACell", summarize_layer="counts")
rna.obs["n_counts"] = np.ravel(rna.X.sum(axis=1))

rna.write(f"{datadir}/TF_regulatory_inference/anndata/lncap_rna_seacells_num_{num}_{date}.h5ad")
SEACell_rna.write(f"{datadir}/TF_regulatory_inference/anndata/lncap_rna_seacell_object_num_{num}_{date}.h5ad")
sc_rna = rna.obs["SEACell"] 
sc_rna.to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_rna_seacells_num_{num}_{date}.csv")

########### Modify based on the soft assignments #########
# If Metacell 1 doesn't fit well (low weight and mismatching sample barcode prefix), and Metacell fits with matching sample and comparable weight, 
# reassign the cell to Metacell 2.

soft_labels, weights = model.get_soft_assignments()

# Create a DataFrame to store the combined information
cell_metacell = pd.DataFrame(index=soft_labels.index)

# Add SEACell labels and weights to the DataFrame interleaved
for i in range(5):
    cell_metacell[f"Metacell_{i+1}_Label"] = soft_labels.iloc[:, i]
    cell_metacell[f"Metacell_{i+1}_Weight"] = weights[:, i]

cell_metacell.to_excel(f"{datadir}/TF_regulatory_inference/outputs/lncap_scRNA_cell_metacell_soft_assignments_num_{num}_{date}.xlsx", index=True)

def update_seacell_labels(cell_metacell, anndata):
    for index, row in cell_metacell.iterrows():
        single_cell_sample_prefix = index.split("_")[0]
        metacell_1_label_prefix = row["Metacell_1_Label"].split("_")[0]
        metacell_2_label_prefix = row["Metacell_2_Label"].split("_")[0]
        
        # IF the sample prefix of the single cell and metacell don't match, check the weight of the second metacell:
        if (single_cell_sample_prefix != metacell_1_label_prefix and
            row["Metacell_1_Weight"] < 0.5 and  # metacell 1 has weight < 0.5
            # the proportion unit between 1st and 2nd metacell is <= 0.2:
            abs(row["Metacell_1_Weight"] - row["Metacell_2_Weight"]) <= 0.2 and
            # change the metacell of the single cell to metacell 2:
            metacell_2_label_prefix == single_cell_sample_prefix):
            anndata.obs["SEACell"] = anndata.obs["SEACell"].astype(str)
            anndata.obs.loc[index, "SEACell"] = row["Metacell_2_Label"]
    
    return anndata

rna_mod = update_seacell_labels(cell_metacell, rna)
rna_mod.write(f"{datadir}/TF_regulatory_inference/anndata/lncap_rna_seacells_mod_num_{num}_{date}.h5ad")

sc_rna_mod = rna_mod.obs["SEACell"]
sc_rna_mod.to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_rna_seacells_mod_num_{num}_{date}.csv")

SEACell_rna_mod = SEACells.core.summarize_by_SEACell(rna_mod, SEACells_label="SEACell", summarize_layer="counts")

SEACell_rna_mod.write(f"{datadir}/TF_regulatory_inference/anndata/lncap_rna_seacell_object_mod_num_{num}_{date}.h5ad")

# re-plot the diagnostics:

# Lower values of compactness suggest more compact/lower variance metacells.
compactness_mod = SEACells.evaluate.compactness(rna_mod, build_kernel_on)
print(compactness_mod.head())
compactness_mod.to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_rna_compactness_mod_num_{num}_{date}.csv")

with PdfPages(f"{datadir}/TF_regulatory_inference/figures/metacell_compactness_mod_num_{num}_{date}.pdf") as pdf:
    plt.figure(figsize=(4, 4))
    sns.boxplot(data=compactness_mod, y="compactness")
    plt.title("Compactness")
    sns.despine()
    pdf.savefig(bbox_inches="tight")
    plt.close()

# Higher values of separation suggest better distinction between metacells.
separation_mod = SEACells.evaluate.separation(rna_mod, "X_pca", nth_nbr=1)
separation_mod.to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_rna_separation_mod_num_{num}_{date}.csv")
with PdfPages(f"{datadir}/TF_regulatory_inference/figures/metacell_separation_mod_num_{num}_{date}.pdf") as pdf:
    plt.figure(figsize=(4, 4))
    sns.boxplot(data=separation_mod, y="separation")
    plt.title("Separation")
    sns.despine()
    pdf.savefig(bbox_inches="tight")
    plt.close()

print("The end")


# save components for the coexpression analysis #
from scipy import io

# first the single-cell data

rna.obs["UMAP_1"] = rna.obsm["X_umap"][:,0]
rna.obs["UMAP_2"] = rna.obsm["X_umap"][:,1]

rna.obs.to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_scRNA_obs_for_coexpression.csv")
rna.var.to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_scRNA_var_for_coexpression.csv")

raw = scipy.sparse.coo_matrix(rna.layers["counts"])
raw = raw.transpose()
from scipy.io import mmwrite
mmwrite(os.path.join(datadir,"/TF_regulatory_inference/utputs/lncap_scRNA_counts_for_coexpression.mtx"), raw)

# save the PCA
pd.DataFrame(rna.obsm["X_pca"]).to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_scRNA_pca_for_coexpression.csv")

SEACell_rna_mod.obs.to_csv(f"{datadir}/TF_regulatory_inference/outputs/lncap_scRNA_seacells_obs_for_coexpression.csv")
rawSc = scipy.sparse.coo_matrix(SEACell_rna_mod.layers["raw"])
rawSc = rawSc.transpose()
mmwrite(f"{datadir}/TF_regulatory_inference/outputs/lncap_scRNA_seacell_counts_for_coexpression.mtx", rawSc)