suppressPackageStartupMessages({
    library(GRaNPA)
    library(qs)
    library(ggplot2)
    library(ComplexHeatmap)
    library(caret)
    library(biomaRt)
    library(dplyr)
    })

set.seed(445566)
indir = "/scratch/svc_td_compbio/users/spsiha/scGRN/glue/process_pruned_GRN/3.0//"
datadir = "/scratch/svc_td_compbio/users/spsiha/scGRN/granpa/"


#### Process the network and the DE genes for GRaNPA #### 

datadir = "/scratch/svc_td_compbio/users/spsiha/scGRN/glue/process_pruned_GRN/3.0/"

grn <- read.csv(paste0(datadir, 'outputs/TF_regulons.csv'))

split_targets <- strsplit(grn$Target.Genes, ",")

# Expand into a long format dataframe
grn <- data.frame(
  TF = rep(grn$Group, sapply(split_targets, length)),
  target = unlist(split_targets),
  stringsAsFactors = FALSE
)
grn_all <- read.csv('/scratch/svc_td_compbio/users/spsiha/scGRN/glue/tf2gene/3.0/regulons_pruned_grn_20240307.csv', stringsAsFactors =T)

grn_all['target']  <- gsub(".*'(.*)'.*", "\\1", grn_all$TargetGene)

# let's put the same weight for TF - target pair
grn_all$weight <- rep(1, nrow(grn_all))
grn$weight <- rep(1, nrow(grn))
head(grn)
grn_all$TargetGene <- NULL

gtf_data <- read.table("/scratch/svc_td_compbio/tools/cellranger_refdata/refdata-gex-GRCh38-2020-A/genes/genes.gtf", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
gtf_data$gene_symbol <- sub(".*gene_name (.*?);.*", "\\1", gtf_data$V9)
gtf_data$ensembl_id <- sub(".*gene_id (.*?);.*", "\\1", gtf_data$V9)
# Keep unique combinations of gene symbols and Ensembl IDs
gene_map <- unique(gtf_data[, c("gene_symbol", "ensembl_id")])


process_grn <- function(grn, gene_map) {
  # Merge by gene symbol and target columns
  grnf <- merge(grn, gene_map, by.x = "target", by.y = "gene_symbol", all.x = TRUE)
  # Check if there are duplicates introduced following the Ensembl ID addition
  duplicates <- grnf %>%
    group_by(target) %>%
    summarize(n_unique_ensembl_ids = n_distinct(ensembl_id)) %>%
    dplyr::filter(n_unique_ensembl_ids > 1)
  
  # Rename column names for GRaNPA
  colnames(grnf)[colnames(grnf) == "ensembl_id"] <- "gene_name"
  colnames(grnf)[colnames(grnf) == "TF"] <- "feature"
  
  # Remove extra columns
  grnf$X <- NULL
  grnf$X0 <- NULL
  
  # Keep unique combinations of column values
  unique_grnf <- unique(grnf[, c("target", "feature", "gene_name", "weight")])
  
  return(unique_grnf)
}

grn_all <- process_grn(grn_all, gene_map)

# Write the results to CSV files
write.csv(grn_all, paste0(datadir, "outputs/TF_regulons_ensembl_grn_all.csv"), quote = FALSE, row.names = FALSE)

DE <- read.csv(paste0(datadir, 'outputs/lncap_scRNA_dataset_de_genes_vs_dmso.csv'))

DEf<- merge(DE, gene_map, by.x = "names", by.y = "gene_symbol", all.x = TRUE)
duplicates <- DEf %>%
  group_by(names) %>%
  summarize(n_unique_ensembl_ids = n_distinct(ensembl_id)) %>%
  dplyr::filter(n_unique_ensembl_ids > 1)

colnames(DEf)[colnames(DEf) == "ensembl_id"] <- "gene_name"
head(DEf)

# remove extra columns
DEf$X = NULL
DEf$X0 = NULL
DEf$names <- NULL
# some of the gene names will be duplicated if there are multiple ensembl ID:s for the same gene,
# let's keep them all but keep the unique combos of column values
unique_DEf <- unique(DEf[, c("comparison", "gene_name","logfoldchanges")])
# rename necessary columns for GRanPa:
colnames(unique_DEf)[colnames(unique_DEf) == "logfoldchanges"] <- "value"
head(unique_DEf)

write.csv(unique_DEf, paste0(datadir, "outputs/lncap_scRNA_dataset_de_genes_for_granpa.csv"), quote=F, row.names=F)


# DEGs dataframe formatted appropriately for GRaNPA: comparison, ENSEMBL ID, log2FC
DE <- read.csv(paste0(indir, "outputs/lncap_scRNA_dataset_de_genes_for_granpa.csv"))
comparisons <- unique(DE$comparison)

# the GRN. We assigned same weight for each TF - gene pair here. The 'target' column is redundant, gene_name carries
# the relevant ENSEMBL ID.
grn_all <- read.csv(paste0(indir, "outputs/TF_regulons_ensembl_grn_all.csv"), stringsAsFactors = TRUE)

grn_all$target <- NULL


# this will be used for filtering the DE genes
logfc_ths = c(0.5)

# modified plotting function with different color scheme and plot title
plot_important_features = function(model, n = 5, group = "") {
  require(dplyr)
  require(cowplot)
  require(ggplot2)

  varImp(model)[1]$importance %>% as.data.frame() %>% arrange(-Overall) -> model
  model$name = rownames(model)
  model = model %>% as.data.frame() %>% dplyr::top_n(Overall, n = n)

  # Define the minimum and maximum colors
  min_color <- "#264653"
  max_color <- "#F4A261"
  model = model %>% arrange(Overall)
  model$name = factor(model$name, levels = model$name)
  # Create the heatmap with the custom color scale
  p1 = model %>% arrange(Overall) %>%
    ggplot(aes(x = name, y = Overall, fill = Overall)) +
    labs(title = paste0("DEGs: ", group))+
    geom_bar(stat = "identity") + coord_flip() +
    theme(plot.title =element_text(size=12))+
    scale_fill_gradientn(colors = c(min_color, max_color), values = c(0, 1)) +  # Set the color scale
    theme_cowplot()

  return(p1)
}

# modified plotting function with different color scheme and plot title
plot_GRaNPA_result <- function(res, ylim = NULL, group = "") {
  require(ggplot2)
  require(cowplot)

  # Create a factor variable for plotting to replace the logical 'random' column
  res$network_type <- factor(res$random, levels = c(FALSE, TRUE), labels = c("Real Network", "Random Network"))
  print(res) # check the rsquared values
  # Define colors for the bars
  colors <- c("#F4A261","#264653")  # Mild colors for distinction

  # Generate the bar plot
  p <- ggplot(res, aes(x = network_type, y = rsq, fill = network_type)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    scale_fill_manual(values = colors) +
    labs(x = "", y = "R-squared (rsq)", title = paste0("DEGs: ", group)) +
    theme_cowplot() +  # Using the cowplot theme
    theme(plot.title=element_text(size=12),legend.position = "none")  # Remove legend since colors are self-explanatory

  # Set y-axis limits if specified
  if (!is.null(ylim)) {
    p <- p + scale_y_continuous(limits = ylim)
  }
  return(p)
}

# Perform GRaNPA analysis for DE set and save TF importance scores
run_granpa_analysis <- function(cond_filtered, comparison, th, grn_name) {
    # Get the appropriate GRN object
    grn <- get(grn_name)
    # Run the GRaNPA analysis
    res <- GRaNPA::GRaNPA_Analysis(gene_metric_data = cond_filtered, GRN_matrix = grn)
    # Define the output files
    result_file <- paste0(datadir, "rds/granpa_result_", grn_name, "_", comparison, "_log2fc_", th, ".qs")
    plot_file <- paste0(datadir, "figures/granpa_result_", grn_name, "_", comparison, "_log2fc_", th, ".pdf")
    tf_importance_file <- paste0(datadir, "figures/granpa_TF_importance_", grn_name, "_", comparison, "_log2fc_", th, ".pdf")
    # Save the results
    qsave(res, result_file)
    res$result['comparison'] = comparison
    # fetch the important features and save to file:
    varImp(res$models$Real)[1]$importance %>% as.data.frame() %>% arrange(-Overall) -> model
    model$name = rownames(model)
    write.csv(model, paste0(datadir, "outputs/granpa_TF_importance_",grn_name,"_",comparison,"_log2fc", th,".csv"))
    # Plot results
    print(paste("Plotting the results for comparison:", comparison, "GRN:", grn_name))
    plot_GRaNPA_result(res$result, group = comparison)
    ggsave(plot_file, width = 4, height = 4)
    # Plot important features
    plot_important_features(res$models$Real, n = 5, group = comparison)
    ggsave(tf_importance_file, width = 4, height = 5)
}


for (comparison in comparisons) {
    cond <- DE[DE$comparison == comparison,]
    cond$comparison <- NULL
        for (th in logfc_ths) {
        # Filter by threshold
        cond_filtered <- cond[abs(cond$value) > th,]
        print(paste("Performing the analysis for comparison:", comparison, "with log2FC threshold:", th))
        # Check if enough genes remain before running analysis
        matching_genes <- sum(cond_filtered$gene_name %in% grn_all$gene_name)
        if (matching_genes < 20) {
            print(paste("Skipping comparison:", comparison, "log2FC:", th, "- Not enough matching genes (", matching_genes, ")"))
            next  # Skip this iteration and move to the next 'th'
        }
        # Call the function to perform the analysis for each GRN
        run_granpa_analysis(cond_filtered, comparison, th, 'grn_all')
    
    }
}

#         rsq ngenes random   comparison   network_type
# 1 0.22583397    183  FALSE resB_vs_dmso   Real Network
# 2 0.02739663    183   TRUE resB_vs_dmso Random Network

#         rsq ngenes random   comparison   network_type
# 1 0.18650906    212  FALSE resA_vs_dmso   Real Network
# 2 0.05610071    212   TRUE resA_vs_dmso Random Network

#         rsq ngenes random     comparison   network_type
# 1 0.31902816    135  FALSE enz48h_vs_dmso   Real Network
# 2 0.06827175    135   TRUE enz48h_vs_dmso Random Network