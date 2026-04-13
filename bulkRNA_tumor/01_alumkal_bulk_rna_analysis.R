suppressPackageStartupMessages({
library(DESeq2)
library(purrr)
library(rtracklayer)
library(survival)
library(survminer)
library(parallel)
library(jsonlite)
library(limma)
library(GSVA)
library(GSEABase)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(SummarizedExperiment)
library(qs)
library(openxlsx)
})


datadir = getwd()

source(paste0(datadir, '/utils.R'))


set.seed(556677)
# Load Ensembl annotation
ann <- getAnnotationTable()
head(ann)

alumkal <- read.table(file.path(datadir, '/data/Alumkal_pnas.1922207117.sd01.csv'), 
                   sep = ";", header=FALSE)


survival_info <- data.frame(read.xlsx(file.path(datadir, '/data/Alumkal_valdata_6-29-2020.xlsx')))
colnames(survival_info) <- gsub("\\s+", "_", colnames(survival_info))
survival_info$os_event <- ifelse(survival_info$os_event == "Yes", 1, 0)

# Convert os_time_mos to numeric
survival_info$os_time_mos <- as.numeric(survival_info$os_time_mos)

# fetch the column names from first row, ignoring the gene_names in the first column
column_names <- as.character(alumkal[1,-1])

alumkal <- alumkal[-1,]
numeric_cols <- colnames(alumkal)[-1]
alumkal[numeric_cols] <- lapply(alumkal[numeric_cols], as.character)
dim(alumkal)

# Replace commas with dots in all numeric columns
alumkal[numeric_cols] <- lapply(alumkal[numeric_cols], function(x) gsub(",", ".", x))
# Convert the columns back to numeric
alumkal[numeric_cols] <- lapply(alumkal[numeric_cols], as.numeric)
# the few duplicates are removed:
# sum(duplicated(alumkal[,1]))
alumkal <- alumkal[!duplicated(alumkal[,1]),]

rownames(alumkal) <- alumkal[,1]
alumkal[,1] <- NULL

colnames(alumkal) <- column_names
alumkal = as.matrix(alumkal)

head(alumkal)

write.table(alumkal, paste0(datadir, '/bulkRNA_tissue/outputs/alumkal_count_matrix.csv'), quote=FALSE, row.names=TRUE, col.names=TRUE)


mat <- alumkal
sample_names <- colnames(alumkal)

# Creating sample_info dataframe
# Define which samples correspond to which group based on the column names
sample_info <- data.frame(
  sample_name = sample_names,
  status = ifelse(1:length(sample_names) <= 7, "non_responder", "responder"),  # Assign status based on index positions
  row.names = sample_names
)
sample_info
# keep only samples assigned as responder or non-responder
sample_info <- sample_info[!is.na(sample_info$status),]
mat <- mat[,sample_info$sample_name]
dim(sample_info)

## Generate gene ranges
rowranges <- ann[ann$gene_name %in% rownames(mat)]
## Filter mithocondrial genes and keep protein coding genes
rowranges <- rowranges[rowranges$gene_type == "protein_coding" & !(seqnames(rowranges) %in% "chrM"), ]

rows_match <- match(rownames(mat), rowranges$gene_name)
rows_to_keep <- !is.na(rows_match)

mat <- mat[rows_to_keep, ]

mat <- log2(mat + 1) 

sample_groups <- factor(c(rep("Non-Responder", 7), rep("Responder", 18)))

gs= f"{datadir}/data/genesets_to_plot.xlsx"

gene_sets <- lapply(gene_sets, function(x) {
  unique(x[!is.na(x)])
})



gene_sets_to_plot <- c('ARS_stress_linked_sensitive', 'SOX4_linked_resistance',
			'SOX4_regulon')

gene_sets_list <- gene_sets_list[names(gene_sets_list) %in% gene_sets_to_plot]


ar_genes <- gene_sets$hallmark_androgen_response
ar_genes <- ar_genes[!is.na(ar_genes)]

gene_sets_list$ARS_stress_linked_sensitive <- setdiff(gene_sets_list$S_stress_linked_sensitive, ar_genes)

# append SOX4 to its regulon
gene_sets_list[['SOX4_regulon']] <- c(
  gene_sets_list[['SOX4_regulon']],
  'SOX4'
)



# Convert to GeneSet objects
gene_sets_list  <- lapply(names(gene_sets_list), function(set_name) {
  GeneSet(unique(gene_sets_list[[set_name]]), setName = set_name)
})

# Build the collection
gene_sets_list  <- GeneSetCollection(gene_sets_list)
mat <- as.matrix(mat)

# Prepare GSVA input
gsvaPar <- gsvaParam(mat, gene_sets_list, kcdf='Gaussian') # Gaussian is the default kernel, Poisson used for integer values
gsva_results <- gsva(gsvaPar)



gsva_results <- gsva_results[, sample_info$sample_name, drop = FALSE]
print('computing the statistics')


# Melt the GSVA results to long format for easy plotting
gsva_long <- reshape2::melt(gsva_results)
colnames(gsva_long) <- c("gene_set", "Sample", "GSVA_Score")
gsva_long['group'] = rep("group", nrow(gsva_long))
# Merge GSVA results with sample information (response status)
gsva_long <- merge(gsva_long, sample_info, by.x = "Sample", by.y = "sample_name")

# Create design matrix for linear model
design <- model.matrix(~ sample_info$status)  # Model with response status
fit <- lmFit(gsva_results, design)  # Fit the linear model
fit <- eBayes(fit)  # Apply empirical Bayes smoothing

diff_enrich <- topTable(fit, coef = "sample_info$statusresponder", number = Inf)
diff_enrich$gene_set <- rownames(diff_enrich)


gsva_long <- merge(gsva_long, 
                   diff_enrich[, c("gene_set", "P.Value")], 
                   by = "gene_set",
                   all.x = TRUE)
 

write.csv(diff_enrich, file = paste0(datadir, '/bulkRNA_tissue/figures/alumkal_limma_p_values.csv'))

# assign 'group' for plotting purposes

gsva_long['group'] = rep('group', nrow(gsva_long))

print('starting to plot')
pdf(paste0(datadir, "/bulkRNA_tissue/figures/alumkal_all_gene_sets_boxplots.pdf"))
for (gene_set_name in unique(gsva_long$gene_set)) { 
  # Subset for the current gene_set
  gsva_long_sub <- gsva_long[gsva_long$gene_set == gene_set_name, ]
  # Create boxplot with individual points
  p <- ggplot(gsva_long_sub, aes(x = group, y = GSVA_Score, fill = status)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 12)
    ) +
    labs(x = "Response Status", y = "GSVA Score", title = gene_set_name) +
    scale_fill_manual(values = c("non_responder" = "#9757C3", "responder" = "#376524")) +
    scale_color_manual(values = c("non_responder" = "#9757C3", "responder" = "#376524")) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    geom_text(
      data = gsva_long_sub,
      aes(label = paste0("p = ", round(P.Value, 3))),	
      x = 1.5,
      y = 1.9,
      size = 3
    )
  print(p)  
}
dev.off()


names(gsva_long)[names(gsva_long) == 'Sample'] <- 'patient_id'
df <- merge(gsva_long, 
                   survival_info, 
                   by = "patient_id",
                   all.x = TRUE)


pdf(paste0(datadir, "/bulkRNA_tissue/figures/Alumkal_all_gene_sets_kaplan_meier.pdf"))
for (gs in as.character(unique(df$gene_set))) {
  df_sub <- df[df$gene_set == gs, ] 
  # Compute median and assign groups
  median_gsva <- median(df_sub$GSVA_Score)
  df_sub$gsva_group <- ifelse(df_sub$GSVA_Score > median_gsva, "High score", "Low score")
  # Survival object
  surv_obj <- Surv(time = df_sub$os_time_mos, event = df_sub$os_event)
  km_fit <- survfit(surv_obj ~ gsva_group, data = df_sub)
  # Kaplan-Meier plot
  km_plot <- ggsurvplot(
    km_fit, 
    data = df_sub, 
    pval = TRUE,
    legend.labs = c("High score", "Low score"),
    palette = c("#E0115F", "#6495ED"),
    risk.table = TRUE,
    title = gs,
    subtitle = subtitle_text
  )
  print(km_plot)   # add to the PDF
}
dev.off()
