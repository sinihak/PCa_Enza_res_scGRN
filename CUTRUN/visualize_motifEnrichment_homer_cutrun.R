suppressPackageStartupMessages({
    library(ggplot2)
})

datadir = getwd()

# read the known results for each peak class
ar = read.delim(paste0(datadir, '/motifs/AR_only_peaks_consensus_filtered_motifs/knownResults.txt'))
nfic = read.delim(paste0(datadir, '/motifs/NFIC_only_peaks_consensus_filtered_motifs/knownResults.txt'))
shared = read.delim(paste0(datadir, '/motifs/Shared_peaks_consensus_filtered_motifs/knownResults.txt'))
head(ar)

# fetch all unique motif occurences
all_motifs = unique(c(ar$Motif.Name, nfic$Motif.Name, shared$Motif.Name))

# list the motif enrichment results
motif_list <- list(AR_only = ar,NFIC_only = nfic,Shared = shared)

# define column names, R messes up the formatting
columns = c("Motif Name", "Consensus", "P-value", "Log P-value", "Benjamini q-value", 
"N of Target Sequences with motif", "Percent of Target Sequences with Motif", "N Background Sequences with Motif",
"Percent of Background Sequences with Motif")

# create empty table, motifs as rows and columns for each peak class
motifs_df =as.data.frame(matrix(nrow=length(all_motifs), ncol=4*length(motif_list)), stringsAsFactors=F)

# assign the motifs as row names
rownames(motifs_df) <- all_motifs

# construct a single motif summary object
for (i in seq_along(motif_list)) {
    df <- motif_list[[i]]
    colnames(df) <- columns
    # remove the "%" mark and ensure numeric
    df[,7] <- as.numeric(gsub("%.*",'', df[,7])) # percentages to numeric
    df[,9] <- as.numeric(gsub("%.*",'', df[,9]))
    df <- df[, c(1,4,5,7,9)]    # keep only the relevant metrics
    matches <- match(rownames(motifs_df), df[,1]) # match the indices across peak classes, missing motifs become NAs
    motifs_df[, (4*(i-1)+1):(4*i)] <- df[matches, 2:5] # fill the master table. For each peak class, insert the values into appropriate columns
    # renaming the columns, each column will be of format peak_class;column_name (e.g. AR_only;Q-value)
    colnames(motifs_df)[(4*(i-1)+1):(4*i)] <-
        paste(
            names(motif_list)[i],
            c("Log P-value", "Q-value", "TargetPct", "BackgroundPct"),
            sep = ";"
        )
}

# final data will be of format motif x peak class data
rownames(motifs_df)=all_motifs


# keep only the motifs of interest for visualization:
moi = unique(read.table(paste0(datadir, "/motifs/cutrun_homer_known_motifs_to_plot.txt"))[,1])
motifs_sub = motifs_df[rownames(motifs_df) %in% moi,]


# this generates a list of long format dataframes for each peak class with fold enrichment and peak class information column
plot_list <- list()

for (i in names(motif_list)) {
    cols <- c(paste0(i, ";Log P-value"), paste0(i, ";Q-value"),paste0(i, ";TargetPct"),paste0(i, ";BackgroundPct"))
    tmp <- motifs_sub[, cols] # tmp contains data only for one peak class at a time
    colnames(tmp) <- c("Log P-value","Q-value","TargetPct","BackgroundPct")
    # calculate the fold enrichment for each motif
    tmp$FoldEnrichment <- tmp$TargetPct / tmp$BackgroundPct
    tmp$Motif <- rownames(motifs_sub) # retain motif names
    tmp$PeakClass <- i # assign peak class as a new column
    plot_list[[i]] <- tmp
}

# convert to dataframe 
plot_data <- do.call(rbind, plot_list)


# qval = (plot_data[,'Q-value'])
# qval[is.na(qval)]=1 #if q-value was missing setting to 1

# if the % of Targets Sequences with Motif is missing, replace with a small value
plot_data$TargetPct[is.na(plot_data$TargetPct)]=0.1

# remove the GEO accession and the /Homer prefix from motif names 
plot_data$Motif <- gsub("\\(GSE[0-9]+\\)|/Homer.*", "", plot_data$Motif)

# remove everything after the first underscore
# plot_data$Motif <- sub("/.*", "", plot_data$Motif)

plot_data$logP <- abs(as.numeric(plot_data$'Log P-value'))

# options(scipen = 999)
# head(plot_data$logP)
# summary(plot_data$logP)
# [1] 2615 2410 2410 2399 2358 2277
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     0.0   285.6  1177.0  1714.2  2407.2 15830.0 
# cap the large logP values (very small P values)

# Colour will represent motif significance (|logP|). Values exceeding 500 were capped to improve visual contrast.
# |logP| = 500 corresponds to P = 10^-500, thus very significant still


# -----------------------------------------------
# Motif enrichment dotplot (Supplementary Fig 6A) 
# -----------------------------------------------

plot_data$logP_plot <- pmin(plot_data$logP, 500)
palette = c("lightgrey", "#FF9933", "#CC0000")

# Sort alphabetically by motif
plot_data_sort <- with(plot_data,  plot_data[order(Motif) , ])

# convert to factor, forcing reversed order (first in alphabet will be on top in plot)
plot_data_sort$Motif <- factor(plot_data_sort$Motif, levels = rev(unique(plot_data_sort$Motif)))

p <- ggplot(plot_data_sort ,aes(x = PeakClass, y = Motif)) +
    geom_point(aes(size = FoldEnrichment,colour = logP_plot))+
    scale_size(range = c(1, 10)) +
    scale_color_gradientn(colors=palette) +
    theme_classic(base_size=16) +
    labs(size = 'Fold Enrichment', color="|Homer Log P|")

pdf(paste0(datadir, "/figures/dotplot_homer_motifs_by_peak_class.pdf"),width = 10, height = 6); p; dev.off()



