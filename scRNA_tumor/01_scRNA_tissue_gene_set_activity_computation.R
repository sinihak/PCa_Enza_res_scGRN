suppressPackageStartupMessages({
    library(qs)
    library(ggplot2)
    library(data.table)
    library(readxl)
    library(AUCell)
    library(Matrix)
    library(Seurat)
    library(GSEABase)
    library(stringr)
    library(msigdbr)
})


datadir = getwd()


#####  create gene sets to test:

genesets = data.frame(read_excel(paste0(datadir, "/data/genesets_to_plot.xlsx")))
genesets_list <- as.list(genesets)
genesets_list[] <- lapply(genesets_list, function(x) x[!is.na(x)]) # remove NA values

#ignore canonical AR target genes and ribosomal genes in ARS- and stress-linked sensitivity:
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark  <- split(hallmark_genesets$gene_symbol, hallmark_genesets$gs_name)
ar_genes <- hallmark[names(hallmark) %in% 'HALLMARK_ANDROGEN_RESPONSE']
genesets_list$'ARS_stress_linked_sensitive' <- setdiff(genesets_list$'ARS_stress_linked_sensitive',ar_genes)

genesets_list[['ARS_stress_linked_sensitive']] <- genesets_list[['ARS_stress_linked_sensitive']][
    !grepl('^(RPL|RPS)', genesets_list[['ARS_stress_linked_sensitive']])
  ]

# remove overlaps between SOX4 regulon and SOX4-linked resistance
genesets_list$'SOX4_linked_resistance_wo_SOX4_regulon' <- setdiff(genesets_list$'SOX4_linked_resistance',network_list$SOX4_regulon)

# We will only use Liu_PTEN_loss_
to_keep = c('SOX4_regulon','SOX4_linked_resistance_wo_SOX4_regulon','SOX4_linked_resistance' 'Liu_PTEN_loss_UP', 'ARS_stress_linked_sensitive','AR_activity')
genesets_list <- genesets_list[names(genesets_list) %in% to_keep]

gs <- list()
for (i in seq_along(genesets_list)) {
    gs[[i]] <- GeneSet(unique(genesets_list[[i]]), setName=names(genesets_list)[i])
}

names(gs) <- names(genesets_list)

runAUCell <- function(raw, gs, dataset_name, geneset_name) {
  # Build rankings
  rankings <- AUCell_buildRankings(raw, plotStats=FALSE)
  # Initialize empty lists for AUC and assignments
  auc <- list()
  auc_scores <- list()
  assign <- list()
  # Run AUCell for each gene set
  for (i in names(gs)) {
    auc[[i]] <- AUCell_run(raw, gs[[i]], aucMaxRank=nrow(rankings)*0.1)
    auc_scores[[i]] <- data.frame(t(getAUC(auc[[i]])))
    assign[[i]] <- AUCell_exploreThresholds(auc[[i]], plotHist=TRUE, assign=TRUE)
  }

  print("saving the cell rankings and rds files")
  # Save the rankings, AUC, and assignment results
  qsave(rankings, paste0(aucdir, "rds/",dataset_name,"_cell_rankings.qs"))
  qsave(auc, paste0(aucdir, "rds/", dataset_name, "_auc.qs"))
  qsave(assign, paste0(aucdir, "rds/", dataset_name,"_aucell_assignment.qs"))


  # Write the results to files
  write.table(auc_scores, file = paste0(aucdir, "outputs/auc_scores_",dataset_name,"_cell_rankings.csv"), sep = ",", col.names = TRUE, row.names = TRUE)
  return(list(auc_scores = auc_scores, rankings = rankings, assign = assign))
}

#######  Lyu et al scRNA-seq data: #######

# fetch the rna count matrix
raw_lyu <- t(readMM(paste0(datadir,'/scRNA_tumor/outputs/lyu_scRNA_geneMatrix.mtx')))
cells_lyu <- read.csv(paste0(datadir,'/scRNA_tumor/outputs/lyu_scRNA_cells.csv'))[,2]
genes_lyu <- read.csv(paste0(datadir,'/scRNA_tumor/outputs/lyu_scRNA_features.csv'))[,2]

colnames(raw_lyu) <- cells_lyu
rownames(raw_lyu) <- genes_lyu

raw_lyu <- as(raw_lyu, "CsparseMatrix")

runAUCell(raw_lyu, gs, dataset_name = 'Lyu_scRNA')

#######  Dong et al scRNA-seq data: #######

# fetch the rna count matrix
raw_dong <- t(readMM(paste0(datadir,'/scRNA_tumor/outputs/dong_scRNA_geneMatrix.mtx')))
cells_dong <- read.csv(paste0(datadir,'/scRNA_tumor/outputs/dong_scRNA_cells.csv'))[,2]
genes_dong <- read.csv(paste0(datadir,'/scRNA_tumor/outputs/dong_scRNA_features.csv'))[,2]

colnames(raw_dong) <- cells_dong
rownames(raw_dong) <- genes_dong

raw_dong <- as(raw_dong, "CsparseMatrix")

runAUCell(raw_dong, gs, dataset_name = 'Dong_scRNA')
