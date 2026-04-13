suppressPackageStartupMessages({
library(TFBSTools)
library(motifmatchr)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(chromVARmotifs)
library(universalmotif)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(patchwork)
library(qs)
library(matrixStats)
library(tidyr)
library(dplyr)

})


datadir <- getwd()


# tar -xvf H14CORE_pwm.tar.gz
# for file in *pwm; do echo "" >> "$file"; done #  adding new line to the end of each pwm file

## Params
cutoff <- 0.0001 ## Motif match p-value cutoff
width <- 7 ## Minimum motif length

pathToPWMs = paste0(datadir, "/data/HOCOMOCO/H14CORE/pwm/")
pwmfiles = list.files(pathToPWMs)

peaks = read.csv(paste0(datadir,"/TF_footprinting/data/SOX4_regulatory_region_of_interest_150kb_peaks.bed"), sep="\t", header=FALSE)
peaks$V6 <- NULL


# select regions of interest
poi <- peaks[,c('V1', 'V2', 'V3')]
colnames(poi) = c('chrom', 'start', 'end')

poi$regulatory_region[poi$chrom == "chr6"] <- "SOX4"

poi$peak <- paste0(poi$chrom,':',poi$start,'-',poi$end)
poi <- poi[poi$regulatory_region == 'SOX4',]

# convert peaks of interest to genomicRanges
poi_ranges <- GRanges(
  seqnames = poi$chr,
  ranges = IRanges(start = poi$start, end = poi$end),
)

# process the pwms
# Convert to motif objects and combine into a PFMatrixList 
empty.list <- list()
for (x in 1:length(pwmfiles)) {
  test1 = universalmotif::read_matrix(paste0(pathToPWMs,pwmfiles[x]), 
                    headers = ">", sep = "", positions = "rows") # universalmotif object
  lel = convert_motifs(test1, class = "TFBSTools-PFMatrix") # convert to PFMatrix 
  empty.list <- c(empty.list, lel) # append to list
}

pfm.list <- do.call(PFMatrixList, empty.list) # create a PFMatrixList from the individual PFMatrix objects
names(pfm.list) <- pwmfiles
motifs <- pfm.list

# get the motif scores for atac peaks
motifscores <- matchMotifs(
  pwms = motifs,
  subject = poi_ranges,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  out = "scores",
  p.cutoff = cutoff,
  w = width
)

# check how many motif matches were obtained:
# print(sprintf("Matches with cutoff=%s and width=%s: %.2f%%", cutoff, width, 100*mean(assay(motifscores,"motifMatches")>0)))

motif_score_matrix <- assay(motifscores, "motifScores")  # Logical matrix of motif matches: peaks x motifs
motif_match_matrix <- assay(motifscores, "motifMatches")  # Logical matrix of motif matches: peaks x motifs

rownames(motif_score_matrix) = poi$peak
rownames(motif_match_matrix) = poi$peak

tf_names <- sub("\\..*$", "", dimnames(motifscores)[[2]])
colnames(motif_score_matrix) <- tf_names
colnames(motif_match_matrix) <- tf_names


motifs_of_interest = c('NFIC','NFIB', 'NFIX','NFIA','FOXA1','ANDR', 'MYC', 'NR3C2', 'NR3C1','CTCF', 'RREB1','ATF4', 'CEBPG', 'SREBF1','MAFK','MAF', 'MAFG','MAFB')

# motif match matrix to bed file 
sub_match_mat <- motif_match_matrix[, colnames(motif_match_matrix) %in% motifs_of_interest,  drop=FALSE]

motif_match_df <- data.frame(sub_match_mat* 1) # converts logical matches to numeric 
# add regulatory region information
motif_match_df$regulatory_region <- poi$regulatory_region

# keep only regulatory region(s) of interest
sox4_matches = motif_match_df[motif_match_df$regulatory_region == 'SOX4',]

df <- sox4_matches
df$peak <- rownames(df)

# parse genomic coordinates (chr:start-end to chrom, start and end columns for BED format)
df <- df %>%
  separate(peak, into = c("chrom", "coords"), sep = ":") %>%
  separate(coords, into = c("start", "end"), sep = "-") %>%
  mutate(
    start = as.numeric(start),
    end   = as.numeric(end)
  )

# convert to long format (peak x TF)
long_df <- df %>%
  pivot_longer(
    cols = -c(chrom, start, end, regulatory_region),
    names_to = "TF",
    values_to = "match"
  )

# keep only motif hits
bed_df <- long_df %>% filter(match == 1)

# create BED-like output
bed_df <- bed_df %>%
  select(chrom, start, end, TF)


# Keep ANDR = AR
andr_bed <- bed_df[grepl("^ANDR", bed_df$TF), ]



write.table(bed_df,file = paste0(datadir,"/TF_footprinting/data/SOX4_region_motif_matches_hocomoco.bed"),  sep = "\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(andr_bed,file = paste0(datadir,"/TF_footprinting/data/SOX4_region_motif_matches_AR_hocomoco.bed"),  sep = "\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)
# ---> andr_bed is used for plotting the genomic track in Fig 4F


### visualize the motif score/match heatmap
# process motif scores
sub_score_mat <- motif_score_matrix[, colnames(motif_score_matrix) %in% motifs_of_interest,  drop=FALSE]
# setdiff(colnames(sub_score_mat), motifs_of_interest)

motif_scores_df <- data.frame(sub_score_mat)
motif_scores_df$regulatory_region <- poi$regulatory_region


regulatory_region <- unique(motif_scores_df$regulatory_region)

# calculate mean motif scores for each region
mean_motif_scores <- sapply(regulatory_region, function(region) {
  subset <- motif_scores_df[motif_scores_df$regulatory_region == region, ]
  subset$regulatory_region <- NULL
  colMeans(subset)
})


# If same TF has multiple matches within the region, calculate sum of their motif scores 
# First fetch the base name of TF (without suffixes)

tf_base <- sub("\\.[0-9]+$", "", rownames(mean_motif_scores))
summed_scores <- as.matrix(tapply(mean_motif_scores[, "SOX4"], tf_base, sum))

# sum the total motif matches for each region
sum_motif_matches <- sapply(regulatory_region, function(region) {
  subset <- motif_match_df[motif_match_df$regulatory_region == region, ]
  subset$regulatory_region <- NULL
  colSums(subset)
})

# Fetch the base name of TF (without suffixes)
tf_base <- sub("\\.[0-9]+$", "", rownames(sum_motif_matches))
# If same TF has multiple matches within the region, calculate their sum 
summed_matches <- as.matrix(tapply(sum_motif_matches[, "SOX4"], tf_base, sum))

# -----------------------------------------------------------------
# SOX4 regulatory region motif match heatmap (Supplementary Fig 5A)
# -----------------------------------------------------------------

p2 <- Heatmap(
  summed_scores,
  col = colorRamp2(c(0, 8), c("white", "darkred")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  name="Motif Score",
  show_row_names = TRUE,
  show_column_names = TRUE,
    # Annotate each tile with sum_motif_matches
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      label = summed_matches[i, j],
      x = x, y = y,
      gp = gpar(fontsize = 10)
    )
  }
)
pdf(paste0(datadir,"/TF_footprinting/plots/LNCaP_scRNA_TFs_of_interest_chromvar_peaks_motif_score_heatmap_SOX4_hocomoco.pdf"), width=4, height=4);draw(p2);dev.off()