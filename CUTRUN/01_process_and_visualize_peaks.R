suppressPackageStartupMessages({
    library(qs)
    library(csaw)
    library(GenomicRanges)
    library(edgeR)
})

datadir = getwd()

temp = list.files(path = paste0(datadir, "/peaks/processed"), pattern="\\.bed$", full.names=T)

# read the peak summits. The replicate peak calls have been concatenated.
peak_summits = lapply(temp, read.delim)
colnames <- c("chr", "start", "end", "peak_name", "score")
peak_summits <- lapply(peak_summits, setNames, colnames)
head(peak_summits[[1]])

peak_summits_gr <- lapply(peak_summits, GRanges)

samples = c("LNCaP_48h_enza_NFIC", "LNCaP_AR", "LNCaP_NFIC", "ResA_NFIC", "ResB_NFIC")
names(peak_summits_gr) <- samples
names(peak_summits) <- samples

# define read parameters
param <- readParam(max.frag = 600, pe="both")

sample_groups <- list(
    LNCaP_48h_enza_NFIC = c("LNCaP_48h_enza_NFIC_1","LNCaP_48h_enza_NFIC_2"),
    LNCaP_AR = c("LNCaP_AR_1", "LNCaP_AR_2"),
    LNCaP_NFIC = c("LNCaP_NFIC_1", "LNCaP_NFIC_2"),
    ResA_NFIC = c("ResA_NFIC_1", "ResA_NFIC_2"),
    ResB_NFIC = c("ResB_NFIC_1", "ResB_NFIC_2")
)

# Compute read counts
peak_counts <- lapply(names(sample_groups), function(sample) {
    message("working on sample ", sample)
    bam_files <- file.path(
        datadir,
        "bams",
        paste0(
            sample_groups[[sample]],
            ".filtered.bl.sorted.bam"
        )
    )
    regionCounts(
        bam_files,
        peak_summits_gr[[sample]],
        param = param
    )
})

names(peak_counts) <- names(sample_groups)

qsave(peak_counts, paste0(datadir, "/peaks/processed/sample_replicate_specific_read_counts.qs"))

# peak_counts <- qread(paste0(datadir, "/peaks/processed/sample_replicate_specific_read_counts.qs"))


# Calculate peak abundances

peak_abundances <- lapply(
    peak_counts,
    function(x) aveLogCPM(asDGEList(x))
)

lapply(peak_counts, function(x) summary(aveLogCPM(asDGEList(x))))
# $LNCaP_48h_enza_NFIC
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.6783  0.7502  1.7299  1.7404  2.6762  7.2467 

# $LNCaP_AR
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.4259  0.3602  1.1327  1.1237  1.8670  6.2770 

# $LNCaP_NFIC
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  -1.250   1.039   1.668   1.742   2.360   6.086 

# $ResA_NFIC
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.7191  1.2891  2.0509  2.0989  2.8635  7.0514 

# $ResB_NFIC
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.4672  1.6787  2.3460  2.3935  3.0660  6.7344 


# --> let's not remove any peaks based on counts yet

# Retain peaks detected in both replicates (count > 0 in both)

peak_summits_filt <- lapply(names(peak_counts), function(s) {
    counts <- assay(peak_counts[[s]])
    keep <- rowSums(counts > 0) == 2
    peak_summits[[s]][keep, ]
})
# assign names
names(peak_summits_filt) <- names(peak_counts)

# Checking the dimensions
lapply(peak_summits, function(x) dim(x))
# > lapply(peak_summits, function(x) dim(x))
# $LNCaP_48h_enza_NFIC
# [1] 229557      5
# $LNCaP_AR
# [1] 218671      5
# $LNCaP_NFIC
# [1] 229055      5
# $ResA_NFIC
# [1] 184302      5
# $ResB_NFIC
# [1] 158528      5

lapply(peak_summits_filt, function(x) dim(x))
# $LNCaP_48h_enza_NFIC
# [1] 213461      5
# $LNCaP_AR
# [1] 214454      5
# $LNCaP_NFIC
# [1] 229019      5
# $ResA_NFIC
# [1] 182406      5
# $ResB_NFIC
# [1] 157865      5


# Fetch NFIC peaks for further processing
nfi_peaks <- peak_summits_filt[c(
    "LNCaP_48h_enza_NFIC",
    "LNCaP_NFIC",
    "ResA_NFIC",
    "ResB_NFIC"
)]

# assign the sample name to a new column
nfi_peaks <- lapply(names(nfi_peaks), function(x) {
    nfi_peaks[[x]]$sample <- x
    nfi_peaks[[x]]
})


names(nfi_peaks) <- c("LNCaP_48h_enza_NFIC","LNCaP_NFIC","ResA_NFIC","ResB_NFIC")

nfi_all <- do.call(rbind, nfi_peaks)

nfi_all_gr <- GRanges(nfi_all)

tail(nfi_all[, c("peak_name", "sample")])

# merge overlapping genomic intervals
nfi_merged <- reduce(nfi_all_gr)
nfi_peaks_gr <- lapply(nfi_peaks, GRanges)



# The NFIC consensus peaks are defined as merged intervals supported
# by at least two independent NFIC datasets (LNCaP, LNCaP-Enza48h, Res-A, RES-B)
support <- sapply(nfi_peaks_gr, function(x) {
    countOverlaps(nfi_merged, x) > 0
})

# find peaks with overlap in at least two samples
# this is now our NFIC consensus set
nfi_consensus <- nfi_merged[rowSums(support) >= 2]

# take AR peaks
# combine AR and NFIC consensus peaks and merge overlaps 
# these are solely derived from LNCaP AR dataset
ar_gr <- GRanges(peak_summits_filt[["LNCaP_AR"]])
# AR consensus peaks:
ar_consensus <- reduce(ar_gr)
all_peaks <- reduce(c(ar_consensus, nfi_consensus))

ar_support <- countOverlaps(all_peaks, ar_consensus) > 0
nfi_support <- countOverlaps(all_peaks, nfi_consensus) > 0

# assign peak class
# AR only: overlaps only AR consensus peaks
# NFIC only: overlaps only NFIC consensus peaks
# Shared: overlaps both AR and NFIC consensus peaks
peak_class <- rep(NA, length(all_peaks)) 
peak_class[ ar_support & !nfi_support] <- "AR_only" 
peak_class[!ar_support & nfi_support] <- "NFIC_only" 
peak_class[ ar_support & nfi_support] <- "Shared" 
all_peaks$peak_class <- peak_class



bams <- unlist(lapply(names(sample_groups), function(sample) {
    bam_files <- file.path(datadir,"bams",paste0(sample_groups[[sample]],".filtered.bl.sorted.bam"))
}))

# calculate peak counts and abundances on candidate consensus peaks
peak_counts_consensus <- regionCounts(bams, all_peaks,param = param)
peak_abundances <- aveLogCPM(asDGEList(peak_counts_consensus))


tapply(peak_abundances, peak_class, summary)


# $AR_only
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.6529 -0.8250 -0.2544 -0.2122  0.2901  3.5070 

# $NFIC_only
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.0204  0.5236  0.9377  1.0523  1.4663  5.0925 

# $Shared
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.3881  1.3398  2.0024  2.1506  2.8374  6.7597 

# assign peak interval as the peak names
rownames(peak_counts_consensus) <- paste0(seqnames(all_peaks), "-", start(all_peaks), "-",end(all_peaks))


sapply(c(-1,-0.75,-0.5,-0.25, 0, 1, 2), function(th) {
sum(peak_abundances > th)
})

# filter based on abundances to obtain the finalish consensus peak set and remove very low-signal peaks
# we'll later filter the shared peaks further
keep <- peak_abundances > -1
peak_counts_consensus <- peak_counts_consensus[keep, ] # only use peaks with mean abundance logCPM > -1

# apply the same filtering to the peak_class object
peak_class_filt <- peak_class[keep]

# add to rowdata
rowData(peak_counts_consensus)$peak_class <- peak_class_filt

counts <- assay(peak_counts_consensus)
colnames(peak_counts_consensus) <- unlist(sample_groups)

qsave(peak_counts_consensus, paste0(datadir, "/rds/ar_nfi_consensus_peak_counts_rse.qs"))

# peak_counts_consensus <- qread(paste0(datadir, "/rds/ar_nfi_consensus_peak_counts_rse.qs"))


# fetch the shared peaks
peaks_shared <- rownames(peak_counts_consensus)[
rowData(peak_counts_consensus)$peak_class == "Shared"
]


samples_reps <- unlist(sample_groups, use.names = FALSE)

# add the colnames 
colnames(peak_counts_consensus) <- samples_reps
colnames(counts) <- samples_reps


AR_cols <- grep("^LNCaP_AR", colnames(counts))
NFIC_cols <- setdiff(seq_len(ncol(counts)), AR_cols)


# some sanity checks to verigy that signal distribution agrees with assigned peak class
# The AR samples should have highest counts on Shared and AR-only peaks, lowest on NFIC-only
tapply(
    rowMeans(counts[, AR_cols]),
    peak_class_filt,
    summary
)

# $AR_only
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    3.00    7.00   10.50   15.83   17.50  526.50 

# $NFIC_only
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.000   4.500   7.500   8.774  11.500  58.000 

# $Shared
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    3.50   20.50   32.00   45.76   54.50 1193.00 


# The NFIC samples should have highest counts on Shared and NFIC-only peaks, lowest on AR-only
tapply(
    rowMeans(counts[, NFIC_cols]),
    peak_class_filt,
    summary
)


# $AR_only
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.125   2.125   3.875   4.697   6.250  52.875 

# $NFIC_only
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.875   9.875  14.000  18.388  21.125 268.500 

# $Shared
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    3.25   16.38   27.62   44.06   52.50  876.75 



# further sanity check, mean signal per peak (in the same cellular background, i.e. control LNCaP parental cells)
ar_means <- rowMeans(counts[, grep("^LNCaP_AR", colnames(counts))])
nfic_means <- rowMeans(counts[, grep("^LNCaP_NFIC", colnames(counts))])

# AR vs. NFIC enrichment
logFC <- log2((ar_means + 1)/(nfic_means + 1)) # pseudocount 1 prevents division by zero
# positive values indicate AR-biased signal, negative NFIC-biased

# plot the ratio 
png(paste0(datadir, "/figures/AR_NFIC_peak_counts_log2fc_boxplot.png"),width = 600,height = 600)

boxplot(
    logFC ~ peak_class_filt,
    ylab = "log2(AR / NFIC)"
)

abline(h = 0, lty = 2, col = "red")
dev.off()


# calculate the distance between peak summits on the peaks

peaks_shared_gr <- all_peaks[all_peaks$peak_class == 'Shared',]

# add peak ID and width
names(peaks_shared_gr) <- paste0(seqnames(peaks_shared_gr), "-", start(peaks_shared_gr), "-",end(peaks_shared_gr))
peaks_shared_gr$peak_id <- paste0(seqnames(peaks_shared_gr), "-", start(peaks_shared_gr), "-",end(peaks_shared_gr))
peaks_shared_gr$length <- width(peaks_shared_gr)

# read the unprocessed summit files (i.e. those without extending by +/- 250bp)
temp2 = list.files(path = paste0(datadir, "/peaks"), pattern="\\.bed$", full.names=T)
names <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(temp2))
names <- sub("_[^_]+$", "", names)
# read the individual (per replicate) peak summits.
peak_summits_raw = lapply(temp2, read.delim)
colnames <- c("chr", "start", "end", "peak_name", "score")
peak_summits_raw <- lapply(peak_summits_raw, setNames, colnames)
names(peak_summits_raw) <- names

# keep only the AR and NFIC samples
peak_summits_raw <- peak_summits_raw[!grepl("2hAb|H3Kme3", names(peak_summits_raw))]
ar_samples <- names(peak_summits_raw[grepl("AR", names(peak_summits_raw))])
nfi_samples <- names(peak_summits_raw[grepl("NFIC", names(peak_summits_raw))])

# convert to genomic ranges
peak_summits_raw_gr <- lapply(peak_summits_raw, GRanges)

# add the peak id
peak_summits_raw_gr <- lapply(peak_summits_raw_gr, function(x){
   x$peak_id <- paste0(seqnames(x), "-", start(x), "-",end(x))
   x
})

# find the AR summit positions overlapping with the shared peaks. Report the start of the summit and the overlapping peak
ar_df <- data.frame()

for (i in ar_samples){
    ar_summits <- peak_summits_raw_gr[[i]]
    overlaps <- findOverlaps(ar_summits, peaks_shared_gr)
    df <- data.frame(
    sample = i,
    summit_pos = start(ar_summits)[queryHits(overlaps)],
    peak_id = names(peaks_shared_gr[subjectHits(overlaps)])
    )
ar_df <- rbind(ar_df, df)
}

colnames(ar_df) <- c('sample', 'summit_start', 'peak_id')

# find the NFIC summit positions overlapping with the shared peaks. Report the start of the summit and the overlapping peak

nfi_df <- data.frame()

for (j in nfi_samples){
    nfi_summits <- peak_summits_raw_gr[[j]]
    overlaps <- findOverlaps(nfi_summits, peaks_shared_gr)
    df <- data.frame(
    sample = j,
    summit_pos = start(nfi_summits)[queryHits(overlaps)],
    peak_id= names(peaks_shared_gr[subjectHits(overlaps)])
    )
nfi_df <- rbind(nfi_df, df)
}

colnames(nfi_df) <- c('sample', 'summit_start', 'peak_id')

# merge by the peak id and calculate the distance between all the AR and NFI summit start positions for the peaks
dist_df <- merge(ar_df, nfi_df, by='peak_id', suffixes=c('_AR', '_NFIC'))
dist_df$distance <- abs(dist_df$summit_start_AR - dist_df$summit_start_NFIC)

# keep only the pairs with distance max. 100 bp
close_pairs <- dist_df[dist_df$distance <= 100,]
head(close_pairs)

peak_summary <- data.frame()
peak_ids <- unique(close_pairs$peak_id)

# create a summary data frame
peak_summary <- do.call(
    rbind,
    lapply(split(close_pairs, close_pairs$peak_id), function(df) { # split close_pairs by the peak_id to handle each peak individually
        if (nrow(df) == 0)
            return(NULL)
        # some distance metrics for each peak
        min_dist <- min(df$distance) # min distance between summit pairs
        mean_dist <- mean(df$distance) # mean distance of summit pairs
        max_dist <- max(df$distance) # max distance between summit pairs
        # report each unique summit pair within the distance of 100 bp 
        pairs <- unique(
            paste(df$sample_AR, df$sample_NFIC, sep = " vs ")
        )
        data.frame(
            peak_id = as.character(df$peak_id[1]),
            pairs = paste(pairs, collapse = "; "),
            stringsAsFactors = FALSE,
            min_distance = min_dist,
            max_distance = max_dist,
            mean_distance = mean_dist
        )
    })
)


# apply filtering to shared peaks
# we will now keep only the shared peaks where the distance between NFIC and AR summit is max. 100 bp

all_peaks_df <- as.data.frame(all_peaks)
all_peaks_df$peak_id <- paste0(seqnames(all_peaks), "-", start(all_peaks), "-",end(all_peaks))

# apply the initial logCP filtering
all_peaks_df_filt <- all_peaks_df[keep, ]

all_peaks_filtered <- all_peaks_df_filt[
    all_peaks_df_filt$peak_class != "Shared" |
    all_peaks_df_filt$peak_id %in% unique(close_pairs$peak_id),
]
# now the filtering based on the summit distances
keep2 <- rownames(peak_counts_consensus) %in% unique(all_peaks_filtered$peak_id)



peak_counts_consensus <- peak_counts_consensus[keep2, ]



#sanity check
all(
    rownames(peak_counts_consensus[keep2, ]) %in%
    all_peaks_filtered$peak_id
    )

# TRUE
qsave(peak_counts_consensus, paste0(datadir, '/rds/ar_nfi_consensus_peak_counts_shared_peaks_filtered_by_summit_dist_max_100bp_rse.qs'))
qsave(all_peaks_filtered, paste0(datadir, '/rds/all_peaks_w_shared_peaks_filtered_by_summit_dist_max_100bp.qs'))
qsave(peak_summary, paste0(datadir, '/rds/summary_shared_peaks_filtered_by_summit_dist_max_100bp.qs'))

all_peaks_filtered <- qread(paste0(datadir, '/rds/all_peaks_w_shared_peaks_filtered_by_summit_dist_max_100bp.qs'))


bed_df <- all_peaks_filtered

# convert to BED coordinates
bed_df$start <- bed_df$start - 1

write.table(bed_df[, c("seqnames", "start", "end", "peak_id", "peak_class")], paste0(datadir, "peaks/all_peaks_consensus_filtered.bed"), sep="\t", quote=F, row.names=F, col.names=F)

peak_classes_unique <- unique(bed_df$peak_class)

for (i in peak_classes_unique){
    df <- bed_df[bed_df$peak_class == i,]
    write.table(df[, c("seqnames", "start", "end", "peak_id")],
    file=paste0(datadir, '/peaks/',i,'_peaks_consensus_filtered.bed'), sep='\t', quote=F, row.names=F, col.names=F)
}

# > mean(close_pairs$distance)
# [1] 39.63237


# -----------------------------------
# Venn plot of peak classes (Fig 5A) 
# -----------------------------------


library(VennDiagram)
set1 <- all_peaks_filtered$peak_id[all_peaks_filtered$peak_class %in% c('AR_only', 'Shared')]
set2 <- all_peaks_filtered$peak_id[all_peaks_filtered$peak_class %in% c('NFIC_only', 'Shared')]



colors = c("#6082B6", "#549D89")
library(VennDiagram)
library(grid)

set1 <- all_peaks_filtered$peak_id[
  all_peaks_filtered$peak_class %in% c("AR_only", "Shared")
]

set2 <- all_peaks_filtered$peak_id[
  all_peaks_filtered$peak_class %in% c("NFIC_only", "Shared")
]


colors <- c("#6082B6", "#549D89")

pdf(file.path(datadir, "figures", "AR_NFIC_peaks_venn.pdf"), height=3, width=3)

venn.plot <- venn.diagram(
  x = list(AR = set1, NFIC = set2),
  fill = colors,
  col=colors,
  cex = 1,
  filename = NULL
)

grid.draw(venn.plot)

dev.off()

