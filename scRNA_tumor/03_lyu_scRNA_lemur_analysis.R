suppressPackageStartupMessages({
library(lemur)
library(SingleCellExperiment)
library(data.table)
library(ggplot2)
library(scales)
library(Matrix)
library(scater)
library(Seurat)
library(harmony)
library(qs)
})

SEED <- 33882
set.seed(SEED)

datadir = getwd()


# fetch the rna count matrix
mtx <- t(readMM(paste0(datadir,"/scRNA_tumor/outputs/lyu_scRNA_geneMatrix.mtx")))
cells <- read.csv(paste0(datadir,"/scRNA_tumor/outputs//lyu_scRNA_cells.csv"))[,2]
genes <- read.csv(paste0(datadir,"/scRNA_tumor/outputs/lyu_scRNA_features.csv"))[,2]
metadata <- read.csv(paste0(datadir,"/scRNA_tumor/outputs/lyu_scRNA_metadata.csv"))
rownames(metadata) = metadata$X
colnames(mtx) <- cells
rownames(mtx) <- genes
seu <- CreateSeuratObject(counts = mtx, meta.data = metadata, assay = "RNA", genome=hg38)
sce_full <- as.SingleCellExperiment(seu)
qsave(sce_full, paste0(datadir,"/scRNA_tumor/rds/scRNA_Lyu_et_al_sce.qs"))

sce_full <- qread(paste0(datadir,"/scRNA_tumor/rds/scRNA_Lyu_et_al_sce.qs"))

groups_to_keep <- c("PC","HSPC")
sub <- sce_full[,sce_full$type %in% groups_to_keep]
sce = sub[,sub$cell_type == "Epithelial cell"]

library(rhdf5)

# Settings
chunk_size <- c(5000, 5000)
tumor_types <- c("HSPC", "PC")
assayNames(sce)

# fit the lemur model
fit <- lemur(
  sce,
  design = ~ type,
  n_embedding = 15
)

# align with harmony
fit <- align_harmony(fit)

qsave(fit, paste0(datadir,"/scRNA_tumor/rds/scRNA_Lyu_et_al_epithelium_lemur_fit_aligned.qs"))

# DE HSPC vs PC (when delta > 0, expression is higher in HSPC)
de <- test_de(
  fit,
  contrast = cond(type = "HSPC") -
             cond(type = "PC")
  )


qsave(de, paste0(datadir, "/scRNA_tumor/rds/scRNA_Lyu_et_al_epithelium_lemur_fit_de_HSPC_vs_PC.qs"))
de_mat <- assay(de, "DE")
h5_file <- paste0(datadir, "/scRNA_tumor/outputs/scRNA_Lyu_et_al_epithelium_lemur_de_HSPC_vs_PC.h5")
h5createFile(h5_file)
h5write(de_mat, h5_file, "de", chunk_size)
write.csv(colnames(de_mat), paste0(datadir, "/scRNA_tumor/outputs/scRNA_Lyu_cells_lemur.csv"), quote=FALSE)
write.csv(rownames(de_mat), paste0(datadir, "/scRNA_tumor/outputs/scRNA_Lyu_genes_lemur.csv"), quote=FALSE)