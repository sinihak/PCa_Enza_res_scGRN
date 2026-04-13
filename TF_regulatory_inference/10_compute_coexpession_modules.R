suppressPackageStartupMessages({library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(qs)
library(igraph)
library(doParallel)
})

registerDoParallel(cores=4)
theme_set(theme_cowplot())
set.seed(555)
num=250



date="20240307"


raw <- Matrix::readMM(paste0(datadir,'/TF_regulatory_inference/outputs/lncap_scRNA_counts_for_coexpression.mtx'))
X_pca <- read.table(paste0(datadir, '/TF_regulatory_inference/outputs/lncap_scRNA_pca_for_coexpression.csv'), sep=',', header=TRUE, row.names=1)
cell_meta <- read.delim(paste0(datadir, '/TF_regulatory_inference/outputs/lncap_scRNA_obs_for_coexpression.csv'), sep=',', header=TRUE, row.names=1)
gene_meta <- read.table(paste0(datadir, '/TF_regulatory_inference/outputs/lncap_scRNA_var_for_coexpression.csv'), sep=',', header=TRUE, row.names=1)
rna_sc = read.csv(paste0(datadir, '/TF_regulatory_inference/outputs/lncap_rna_seacells_integrated_num_',num,'_',date,'.csv'))

### set up the SEACell object:
raw_sea <- Matrix::readMM(paste0(datadir,'/TF_regulatory_inference/outputs/lncap_scRNA_seacell_counts_for_coexpression.mtx'))
sea_meta <- read.delim(paste0(datadir, '/TF_regulatory_inference/outputs/lncap_scRNA_seacells_obs_for_coexpression.csv'), sep=',', header=TRUE, row.names=1)

colnames(raw_sea) <- rownames(sea_meta)
rownames(raw_sea) <- rownames(gene_meta)
dim(raw_sea)

mad_values <- apply(raw_sea, 1, function(x) mad(x, constant = 1))
zero_mad_genes <- names(mad_values)[mad_values == 0]


raw_sea = raw_sea[!rownames(raw_sea) %in% zero_mad_genes,]

sea <- CreateSeuratObject(counts = raw_sea, assay = "RNA", project = "SEACells", min.features = 0, min.cells = 0)

# add the modality information to be used later for expression mtx setup
# the SEACells were computed across the entire dataset, so we will use an attribute that covers the entire dataset
sea[['modality']] = 'RNA'

# get the umap from cell_meta:
umap <- cell_meta[,c('UMAP_1', 'UMAP_2')]

# set the rownames and colnames for the expression matrix:
# for Seurat, rows of X are genes, cols of X are cells
colnames(raw) <- rownames(cell_meta)
rownames(raw) <- rownames(gene_meta)
rownames(X_pca) <- rownames(cell_meta)
rownames(umap) <- rownames(cell_meta)

#  Convert rownames of sparse matrix to logical index for subsetting
gene_to_remove <- rownames(raw) %in% zero_mad_genes
# ignore the zero MAD genes
raw <- raw[!gene_to_remove,]

gene_meta <- gene_meta[!rownames(gene_meta) %in% zero_mad_genes,] 

# rename X_pca columns to Seurat's format
for (i in 1:50) {
  colnames(X_pca)[i] <- paste0("PC_", i)
}

rna <- CreateSeuratObject(counts = raw,  meta.data = cell_meta,  assay = "RNA",  project = "SEACells", min.features = 0,  min.cells = 0)

# Create the var.features info where FALSE becomes NA
rna@assays$RNA@meta.data$var.features <- ifelse(gene_meta$highly_variable == "True", rownames(gene_meta), NA)

head(rna@assays)

hvgs <- rownames(gene_meta)[gene_meta$highly_variable == "True"]

# get the metacell object
metacell_obj <- GetMetacellObject(rna,wgcna_name = 'SEACells$wgcna_metacell_obj')

#assign the HVGs to the seurat object (although full GEX landscape is used for the module construction)
hvgs <- rownames(gene_meta)[gene_meta$highly_variable == "True"]
VariableFeatures(rna) <- hvgs

# set PCA reduction
rna@reductions$pca <- CreateDimReducObject(embeddings = as.matrix(X_pca), key="PC",  assay="RNA")
# set UMAP reduction
rna@reductions$umap <- CreateDimReducObject(embeddings = as.matrix(umap), key="UMAP", assay="RNA")

# normalize expression matrix
rna <- NormalizeData(rna)

# set up hdWGCNA experiment
rna <- SetupForWGCNA(rna,
    wgcna_name = 'SEACells$wgcna'
)

length(rna@misc$SEACells$wgcna)

# add the seacells dataset
rna <- SetMetacellObject(rna, sea)
rna <- NormalizeMetacells(rna)

# setup expression matrix
rna <- SetDatExpr(
  rna,
  group.by='modality',
  group_name='RNA', # This same column should have also been used for metacell computation (or more like, the same groups should exist in the metacell object)
  use_metacells=TRUE,
  assay='RNA', 
  layer = "data"
)

# expr_matrix <- GetAssayData(metacell_obj, assay = "RNA", layer = "data")[zero_mad_genes, ]
# apply(expr_matrix, 1, function(x) range(x))  # Check min & max expression

qsave(rna, paste0(datadir, '/TF_regulatory_inference/rds/lncap_scRNA_seurat_object_for_coexpression_analysis.qs'))
rna <- qread(paste0(datadir, '/TF_regulatory_inference/rds/lncap_scRNA_seurat_object_for_coexpression_analysis.qs'))


# test soft power threshold
rna <- TestSoftPowers(rna, networkType = 'signed')

# plot the results:
plot_list <- PlotSoftPowers(rna)
# assemble with patchwork
pdf(paste0(datadir, '/TF_regulatory_inference/figures/lncap_scRNA_coexpression_soft_powers.pdf'))
wrap_plots(plot_list, ncol=2)
dev.off()

# compute the co-expression network
rna <- ConstructNetwork(rna,overwrite_tom = TRUE)

# compute module eigengenes and eigengene-based connectivity
rna <- ModuleEigengenes(rna)
rna <- ModuleConnectivity(rna)

# rename modules
rna <- ResetModuleNames(
  rna,
  new_name = 'sc-M',
  wgcna_name='SEACells$wgcna'
)

qsave(rna, paste0(datadir, '/TF_regulatory_inference/rds/lncap_scRNA_coexpression_analysis.qs'))
# rna <- qread(paste0(datadir, 'rds/lncap_scRNA_coexpression_analysis.qs'))
modules <- GetModules(rna, wgcna_name='SEACells$wgcna')
write.csv(modules, paste0(datadir, '/TF_regulatory_inference/outputs/coexpression_modules.csv'),quote=FALSE)
# fetch harmonized module eigengenes 
MEs<- GetMEs(rna, TRUE,wgcna_name='SEACells$wgcna')

write.csv(MEs, paste0(datadir, '/TF_regulatory_inference/outputs/lncap_scRNA_module_eigengenes.csv'),quote=FALSE, row.names=TRUE)

# get WGCNA network and module data
net <- GetNetworkData(rna, wgcna_name="SEACells$wgcna")

module_genes <- modules$gene_name
module_colors <- modules$color
names(module_colors) <- modules$gene_name

pdf(paste0(datadir, "/TF_regulatory_inference/figures/lncap_scRNA_coexpression_dendro.pdf"))
PlotDendrogram(rna)
dev.off()

pdf(paste0(datadir, "/TF_regulatory_inference/figures/lncap_scRNA_eigengene-based_connectivity_rank.pdf"))
PlotKMEs(rna, ncol=2)
dev.off()