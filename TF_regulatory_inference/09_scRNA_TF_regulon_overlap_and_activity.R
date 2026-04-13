suppressPackageStartupMessages({
    library(qs)
    library(ggplot2)
    library(data.table)
    library(AUCell)
    library(readxl)
    library(stringr)
    library(GSEABase)
    library(Matrix)
})

set.seed(321)
datadir = getwd()

date="20240307"
GRNqval='0.1'

# gene sets
df <- read.csv(paste0(datadir, "/TF_regulatory_inference/outputs/regulons_pruned_grn_",date,".csv"))

# extract gene information inside parentheses
df$genes <- str_extract(df[,2], "'(.*?)'")
# remove quotes
df$genes <- gsub("'", "", df$genes)
# create a list of genes for each TF
regulons <- split(df$genes, df$TF)
# remove empty strings
regulons <- lapply(regulons, function(x) x[grep("\\S", x)])
head(regulons)
cat("This list has", length(names(regulons)), "TFs.\n")

cat("This list has", length(names(regulons)), "TFs.\n")

qsave(regulons,paste0(datadir, "/TF_regulatory_inference/outputs/regulons_pruned_grn_as_list_",date,".qs"))

# Calculate the Jaccard index between the regulons
# Supplementary Fig 1D

jacc.mtx <- matrix(NA, nrow = nRegulons, ncol = nRegulons,dimnames = list(names(regulons), names(regulons)))

# Function to compute Jaccard index between two TF regulons
ComputeJaccardIndex <- function(regulon1, regulon2) {
  set1 <- unique(regulon1)
  set2 <- unique(regulon2)
  # Compute Jaccard index
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  # Check if the union is 0 (to avoid division by zero)
  if (union == 0) {
    return(NA)  # Return NA if union is 0
  }
  jaccard_index <- round(intersection / union,3)
  return(jaccard_index)
}

for (i in 1:nRegulons) {
  for (j in 1:nRegulons) {    
    # Compute Jaccard index between regulons i and j
    jacc.idx <- ComputeJaccardIndex(regulons[[i]], regulons[[j]])
    jacc.mtx[i, j] <- jacc.idx
  }
}


# Update row and column names in jacc.mtx
dimnames(jacc.mtx) <- list(regulon_names_with_counts, regulon_names_with_counts)
jacc.mtx

hm1 <- Heatmap(jacc.mtx, name = "Jaccard index", show_row_names = TRUE, show_column_names = TRUE,
                cluster_columns = TRUE, cluster_rows = TRUE, 
                col = colorRamp2(c(0,0.5,1), c("white","#004d99", "#000d1a")),
                row_names_gp = gpar(fontsize = 9),
                column_names_gp = gpar(fontsize = 9),
                row_dend_side = "left", column_dend_side = "top",
                column_title = NULL, row_title = NULL,
                width = unit(9, "cm"),
                height = unit(9, "cm"),
                heatmap_legend_param = list(title = "JI"))


pdf(paste0(datadir, "/TF_regulatory_inference/figures/lncap_heatmap_regulon_jaccard_num_",num,"_",date,".pdf"));draw(hm1);dev.off() 


# Now proceed to calculate AUCell for regulons

gs <- list()
for (i in names(regulons)){
    gs[[i]] <- GeneSet(unique(regulons[[i]]), setName=paste0(i, "_",length(unique(regulons[[i]])),"g"))
}

qsave(gs, paste0(datadir, "/TF_regulatory_inference/rds/TF_regulon_geneset_for_AUCell_",date,".qs"))
# gs <- qread(paste0(datadir, "(rds/TF_regulon_geneset_for_AUCell_",date,".qs"))

# fetch the rna count matrix
raw <- t(readMM(paste0(datadir,"/TF_regulatory_inference/outputs/lncap_scRNA_raw_mtx_",date,".mtx")))
cells <- read.csv(paste0(datadir,"/TF_regulatory_inference/outputs/lncap_scRNA_cells_",date,".csv"))[,2]
genes <- read.csv(paste0(datadir,"/TF_regulatory_inference/outputs/lncap_scRNA_features_",date,".csv"))[,2]


colnames(raw) <- cells
rownames(raw) <- genes

print(class(raw))

raw <- as(raw, "CsparseMatrix")
print(class(raw))

runAUCell <- function(raw, gs, dataset_name, geneset_name) {
  # Build rankings
  rankings <- AUCell_buildRankings(raw, plotStats=TRUE)

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

  # Save the rankings, AUC, and assignment results
  qsave(rankings, paste0(datadir, "/TF_regulatory_inference/rds/",dataset_name,"_",geneset_name,"_cell_rankings_", date, ".qs"))
  qsave(auc, paste0(datadir, "/TF_regulatory_inference/rds/", dataset_name, "_", geneset_name, "_auc_", date, ".qs"))
  qsave(assign, paste0(datadir, "/TF_regulatory_inference/rds/", dataset_name, "_", geneset_name, "_aucell_assignment_", date, ".qs"))

  # Write the results to files
  write.table(auc_scores, file = paste0(datadir, "/TF_regulatory_inference/outputs/auc_scores_",dataset_name,"_",geneset_name,"_", date, ".csv"), sep = ",", col.names = TRUE, row.names = TRUE)
#   return(list(auc_scores = auc_scores, rankings = rankings, assign = assign))
}


runAUCell(raw, gs, dataset_name = "lncap_scRNA",geneset_name="TF_regulons")

# Remove "self-loops", i.e. the occurences of the TF regulating itself

regulons_filtered  <- lapply(names(regulons), function(tf) {
  regulons[[tf]][regulons[[tf]] != tf]
})

names(regulons_filtered) <- names(regulons)

gs <- list()
for (i in names(regulons_filtered)){
    gs[[i]] <- GeneSet(unique(regulons_filtered[[i]]), setName=paste0(i, "_",length(unique(regulons_filtered[[i]])),"g"))
}

runAUCell(raw, gs, dataset_name = "lncap_scRNA",geneset_name="TF_regulons_wo_self_loops")

# Some extra genesets
# First hallmark genes. Overlaps are removed for correlation analysis with SOX4 regulon
gs_hm <- qread(paste0(datadir, "/TF_regulatory_inference/rds/msigdb_hallmark_for_AUCell.qs"))
gs_hm_filt <- qread(paste0(datadir, "/TF_regulatory_inference/rds/msigdb_hallmark_filtered_from_overlaps_w_sox4_regulon_for_AUCell.qs"))

runAUCell(raw, gs_hm, dataset_name = "lncap_scRNA",geneset_name="msigdb_hallmark")
runAUCell(raw, gs_hm_filt, dataset_name = "lncap_scRNA",geneset_name="msigdb_hallmark_overlaps_w_sox4_regulon_filtered")

# XBP1 target genes
genesets = data.frame(read_excel(paste0(datadir, "/data/genesets_to_plot.xlsx"))
genesets_list <- as.list(genesets)
genesets_list[] <- lapply(genesets_list, function(x) x[!is.na(x)])
XBP1_targets <- list(GeneSet(genesets_list[["XBP1s_Activates_Chaperone_Genes_R_HSA_381038"]], setName="XBP1s_Activates_Chaperone_Genes_R_HSA_381038"))
names(XBP1_targets) <- c("XBP1s_Activates_Chaperone_Genes_R_HSA_381038")

runAUCell(raw, XBP1_targets, dataset_name = "lncap_scRNA",geneset_name="XBP1_targets")
