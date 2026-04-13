
suppressPackageStartupMessages({
library(Signac)
library(Seurat)
library(chromVAR)
library(TFBSTools)
library(JASPAR2022)
library(motifmatchr)
library(chromVARmotifs)
library(ComplexHeatmap)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(SummarizedExperiment)
library(dplyr)
library(data.table)
library(patchwork)
library(rtracklayer)
library(ggplot2)
library(qs)
})

set.seed(456)

datadir = getwd()
date='20240307'

motif <- qread(paste0(datadir,'/TF_footprinting/rds/lncap_atac_motif_object_',date,'.qs'))
atac <- qread(paste0(datadir,'/TF_footprinting/rds/lncap_atac_motif_activities_',date,'.qs'))


# compare the resistant to DMSO
atac$condition <- ifelse(
  grepl("res", atac$dataset), "resistant",
  ifelse(
    atac$dataset == "dmso", "parental",
    ifelse(
      atac$dataset == "enz48h", "parental_enza",
      "unknown"
    )
  )
)

atac[['chromvar']]

DefaultAssay(atac) <- 'chromvar'
Idents(atac) <- 'condition'

diff_activity <- FindMarkers(
  object = atac,
  ident.1 = 'resistant',
  ident.2 = 'parental',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)


motif_names <- ConvertMotifID(object = motif, id = rownames(diff_activity))
diff_activity$motif_id <- rownames(diff_activity)


rownames(diff_activity) <- motif_names
qsave(diff_activity, paste0(datadir,'/TF_footprinting/rds/lncap_scATAC_diff_activity_resistant_vs_parental.qs'))

diff_activity <- qread(paste0(datadir,'/TF_footprinting/rds/lncap_scATAC_diff_activity_resistant_vs_parental.qs'))

# # Compute safe -log10
# the values smaller than 1e-300 are capped at 1e-300, after which the base-10 log is taken
# avoids the inf values
diff_activity$log10_padj <- -log10(pmax(diff_activity$p_val_adj, 1e-300))
diff_activity_sig <- diff_activity[diff_activity$p_val_adj < 0.05,]
write.table(diff_activity_sig, paste0(datadir, '/TF_footprinting/outputs/LNCaP_scATAC_res_vs_dmso_diff_act_motifs.tsv'),quote=FALSE)


# top down and upregulated motifs based on 95th and 5th percentiles
diff_activity_up <- diff_activity_sig[diff_activity_sig$avg_diff > 0, ]
diff_activity_down <- diff_activity_sig[diff_activity_sig$avg_diff < 0, ]


q5 <- quantile(diff_activity_sig$avg_diff, probs = 0.05, na.rm = TRUE)
q95 <- quantile(diff_activity_sig$avg_diff, probs = 0.95, na.rm = TRUE)
top_up <- diff_activity_up[diff_activity_up$avg_diff > q95,]
top_down <- diff_activity_down[diff_activity_down$avg_diff < q5,]

down <- top_down[order(top_down$avg_diff, decreasing=FALSE), ]
rownames(down)

up <- top_up[order(top_up$avg_diff, decreasing=TRUE), ]
# up[up$log10_padj == 300.0000,]
# rownames(up)

# diff_activity[rownames(diff_activity) == 'FOXA1',]

# Combine
top_motifs <- rbind(top_up, top_down)
rownames(top_motifs)
# # Add motif names as a column
diff_activity$motif <- rownames(diff_activity)
top_motifs$motif <- rownames(top_motifs)
top_motifs
library(ggplot2)
library(ggrepel)


# --------------------------------------------------
# Differential motif activity volcano plot (Fig. 4B) 
# --------------------------------------------------

# let's visualize the very top values to avoid making a very large plot
diff_activity<- diff_activity[diff_activity$log10_padj >= 175, ]
top_motifs <- top_motifs[top_motifs$log10_padj >= 175, ]

p4 <- ggplot(diff_activity, aes(x = avg_diff, y = log10_padj)) +
  geom_point(alpha = 0.5, color = "grey70", size = 3) +  # background points
  geom_point(data = top_motifs, aes(x = avg_diff, y = log10_padj, color = motif), size = 4) +
  xlim(-4,4) +
  geom_text_repel(data = top_motifs,
                  aes(x = avg_diff, y = log10_padj, label = motif),
                  size = 5,
                  max.overlaps = 50,
                  box.padding = 0.9, # extra space around labels
                  point.padding = 0.6) +
  theme_minimal() +
  labs(
    x = "Average Difference",
    y = expression(-log[10]("Adjusted p-value"))
  ) +
  theme(legend.position = "none")  # Hide legend if too cluttered
# Save to PDF
pdf(paste0(datadir, 'figures/lncap_scATAC_resistant_vs_parental_diff_active_motifs.pdf'), width = 7, height = 8)
print(p4)
dev.off()
