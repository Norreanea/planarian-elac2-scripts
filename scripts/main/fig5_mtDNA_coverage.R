############################################################
# Figure 5 — mtDNA coverage with gene annotation
# Source: GitHub repo
#   - mtDNA_anno_corrected.gtf
#   - folder mtDNA_coverage/ with *coverage.txt files
############################################################

## Minimal deps (install if needed)
if (!requireNamespace("dplyr",           quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr",           quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("ggplot2",         quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("viridis",         quietly = TRUE)) install.packages("viridis")
if (!requireNamespace("hrbrthemes",      quietly = TRUE)) install.packages("hrbrthemes")
if (!requireNamespace("GenomicRanges",   quietly = TRUE)) install.packages("GenomicRanges")
if (!requireNamespace("IRanges",         quietly = TRUE)) install.packages("IRanges")
if (!requireNamespace("GenomicFeatures", quietly = TRUE)) install.packages("GenomicFeatures")
if (!requireNamespace("Gviz",            quietly = TRUE)) install.packages("Gviz")
if (!requireNamespace("rtracklayer",     quietly = TRUE)) install.packages("rtracklayer")

## ------------------------------------------------------------------
## GitHub locations
## ------------------------------------------------------------------
base_raw   <- "https://raw.githubusercontent.com/Norreanea/planarian-elac2-scripts/main/data"
gtf_url    <- file.path(base_raw, "mtDNA_anno_corrected.gtf")
cov_base   <- file.path(base_raw, "mtDNA_coverage")

## ------------------------------------------------------------------
## Declare the coverage files that exist in mtDNA_coverage/
## (edit this vector if your repo gets new files)
## ------------------------------------------------------------------
lanes3 <- c("L13","L23","L33")  # dpa3
lanes5 <- c("L15","L25","L35")  # dpa5
batches <- c("batch1","batch2")
groups  <- c("Elac","WT","GFP")

# Build filenames like "Elac_L13_batch1_coverage.txt", etc.
make_names <- function(group, lanes) {
  as.vector(outer(
    paste0(group, "_", lanes, "_"),
    paste0(batches, "_coverage.txt"),
    FUN = "paste0"
  ))
}

cov_list <- c(
  make_names("Elac", lanes3), make_names("Elac", lanes5),
  make_names("WT",   lanes3), make_names("WT",   lanes5),
  make_names("GFP",  lanes3), make_names("GFP",  lanes5)
)

## Labels without “_coverage.txt”
strip_cov <- function(x) sub("_coverage\\.txt$", "", x)
new_name  <- strip_cov(cov_list)

## ------------------------------------------------------------------
## Read all coverage .txt from GitHub
## Each file has columns: seqnames start end score
## ------------------------------------------------------------------
all_cov_df <- data.frame()
for (j in seq_along(cov_list)) {
  url_j <- file.path(cov_base, cov_list[j])
  tmp   <- utils::read.table(url_j, header = FALSE,
                             col.names = c("seqnames","start","end","score"))
  tmp$Type <- new_name[j]
  all_cov_df <- rbind(all_cov_df, tmp)
}


## ------------------------------------------------------------------
## Collapse batches → ELAC/WT/GFP at dpa3 and dpa5
##   dpa3 → *_L13, *_L23, *_L33 ; dpa5 → *_L15, *_L25, *_L35
## ------------------------------------------------------------------
cov_wide <- tidyr::pivot_wider(
  all_cov_df,
  names_from = "Type", values_from = "score"
)

# Helper to row-sum columns matching a regex
sum_by <- function(df, pattern) {
  mat <- dplyr::select(df, dplyr::matches(pattern))
  if (ncol(mat) == 0) return(rep(0, nrow(df)))
  rowSums(as.data.frame(mat), na.rm = TRUE)
}

cov_wide <- cov_wide |>
  dplyr::mutate(
    ELAC3 = sum_by(dplyr::cur_data_all(), "^Elac_L(13|23|33)_batch"),
    ELAC5 = sum_by(dplyr::cur_data_all(), "^Elac_L(15|25|35)_batch"),
    WT3   = sum_by(dplyr::cur_data_all(), "^WT_L(13|23|33)_batch"),
    WT5   = sum_by(dplyr::cur_data_all(), "^WT_L(15|25|35)_batch"),
    GFP3  = sum_by(dplyr::cur_data_all(), "^GFP_L(13|23|33)_batch"),
    GFP5  = sum_by(dplyr::cur_data_all(), "^GFP_L(15|25|35)_batch")
  ) |>
  dplyr::select(.data$seqnames, .data$start, .data$end, ELAC3, ELAC5, WT3, WT5, GFP3, GFP5)

## ------------------------------------------------------------------
## Optional quick plot — coverage by group (sanity check)
## ------------------------------------------------------------------
cov_long <- tidyr::pivot_longer(
  cov_wide,
  cols = c("ELAC3","ELAC5","WT3","WT5","GFP3","GFP5"),
  names_to = "group", values_to = "score"
) |>
  dplyr::mutate(
    dpa  = ifelse(endsWith(.data$group, "3"), "dpa 3", "dpa 5"),
    gene = dplyr::case_when(
      startsWith(.data$group, "ELAC") ~ "ELAC2 KD",
      startsWith(.data$group, "WT")   ~ "WT",
      startsWith(.data$group, "GFP")  ~ "GFP mock",
      TRUE ~ "Other"
    )
  )

p_cov <- cov_long |>
  dplyr::group_by(.data$seqnames, .data$start, .data$dpa, .data$gene) |>
  dplyr::summarise(score = sum(.data$score, na.rm = TRUE), .groups = "drop") |>
  ggplot2::ggplot(ggplot2::aes(x = start, y = score, color = gene)) +
  ggplot2::geom_line() +
  viridis::scale_color_viridis(discrete = TRUE) +
  hrbrthemes::theme_ipsum() +
  ggplot2::labs(title = "Normalized coverage", y = "Normalized coverage", x = NULL) +
  ggplot2::facet_grid(~ dpa)
print(p_cov)

## ------------------------------------------------------------------
## Build GRanges for dpa3 and dpa5 (ELAC, WT, GFP columns)
## ------------------------------------------------------------------
dpa3_df <- cov_wide |>
  dplyr::transmute(
    seqnames  = .data$seqnames,
    start     = as.integer(.data$start),
    end       = as.integer(.data$start) + 1L,
    ELAC_dpa3 = .data$ELAC3,
    WT_dpa3   = .data$WT3,
    GFP_dpa3  = .data$GFP3
  ) |>
  dplyr::arrange(.data$start)

dpa5_df <- cov_wide |>
  dplyr::transmute(
    seqnames  = .data$seqnames,
    start     = as.integer(.data$start),
    end       = as.integer(.data$start) + 1L,
    ELAC_dpa5 = .data$ELAC5,
    WT_dpa5   = .data$WT5,
    GFP_dpa5  = .data$GFP5
  ) |>
  dplyr::arrange(.data$start)

dpa3_gr <- GenomicRanges::makeGRangesFromDataFrame(dpa3_df, keep.extra.columns = TRUE)
dpa5_gr <- GenomicRanges::makeGRangesFromDataFrame(dpa5_df, keep.extra.columns = TRUE)

## ------------------------------------------------------------------
## Import GTF directly from GitHub and make TxDb + GeneRegionTrack
## ------------------------------------------------------------------
gtf_gr  <- rtracklayer::import(gtf_url)  # mtDNA_anno_corrected.gtf
txdb    <- GenomicFeatures::makeTxDbFromGRanges(gtf_gr)
geneTrack <- Gviz::GeneRegionTrack(txdb, options(ucscChromosomeNames = FALSE))

# Put friendly names on the gene track
g2s <- S4Vectors::mcols(gtf_gr)[, c("gene_id","gene_name")]
g2s <- unique(as.data.frame(g2s))
rownames(g2s) <- g2s$gene_id
S4Vectors::mcols(geneTrack)$symbol <- g2s[S4Vectors::mcols(geneTrack)$gene, "gene_name"]

## ------------------------------------------------------------------
## Gviz tracks: axis + coverage (dpa3/dpa5) + gene models
## ------------------------------------------------------------------
ax     <- Gviz::GenomeAxisTrack()
accDT3 <- Gviz::DataTrack(dpa3_gr, type = "a",
                          groups = c("ELAC_dpa3","WT_dpa3","GFP_dpa3"),
                          name = "Normalized coverage")
accDT5 <- Gviz::DataTrack(dpa5_gr, type = "a",
                          groups = c("ELAC_dpa5","WT_dpa5","GFP_dpa5"),
                          name = "Normalized coverage")

## ------------------------------------------------------------------
## Regions to show (same spirit as before)
## ------------------------------------------------------------------
reg_starts <- c(2350, 4400, 5600, 6250, 7280, 8300, 8750, 11470, 13505, 14450)
reg_ends   <- c(2500, 4700, 5700, 6450, 7330, 8400, 8850, 12200, 13800, 14600)

## Region 1
Gviz::plotTracks(
  list(ax, accDT3, accDT5, geneTrack),
  from = reg_starts[1], to = reg_ends[1],
  showId = TRUE, geneSymbol = TRUE,
  fontcolor.feature = "darkblue",
  name = "Normalized coverage",
  featureAnnotation = "gene_name",
  just.group = "left",
  CDS = "darkred", ncRNA = "darkgreen",
  showFeatureId=TRUE,
  main="Region 1",cex.main = 3
)

## Region 2
Gviz::plotTracks(
  list(ax, accDT3, accDT5, geneTrack),
  from = reg_starts[2], to = reg_ends[2],
  showId = TRUE, geneSymbol = TRUE,
  fontcolor.feature = "darkblue",
  name = "Normalized coverage",
  featureAnnotation = "gene_name",
  just.group = "left",
  CDS = "darkred", ncRNA = "darkgreen",
  showFeatureId=TRUE,
  main="Region 2",cex.main = 3
)

## Region 3
Gviz::plotTracks(
  list(ax, accDT3, accDT5, geneTrack),
  from = reg_starts[3], to = reg_ends[3],
  showId = TRUE, geneSymbol = TRUE,
  fontcolor.feature = "darkblue",
  name = "Normalized coverage",
  featureAnnotation = "gene_name",
  just.group = "left",
  CDS = "darkred", ncRNA = "darkgreen",
  showFeatureId = TRUE,
  main = "Region 3", cex.main = 3
)

## Region 4
Gviz::plotTracks(
  list(ax, accDT3, accDT5, geneTrack),
  from = reg_starts[4], to = reg_ends[4],
  showId = TRUE, geneSymbol = TRUE,
  fontcolor.feature = "darkblue",
  name = "Normalized coverage",
  featureAnnotation = "gene_name",
  just.group = "left",
  CDS = "darkred", ncRNA = "darkgreen",
  showFeatureId = TRUE,
  main = "Region 4", cex.main = 3
)

## Region 5
Gviz::plotTracks(
  list(ax, accDT3, accDT5, geneTrack),
  from = reg_starts[5], to = reg_ends[5],
  showId = TRUE, geneSymbol = TRUE,
  fontcolor.feature = "darkblue",
  name = "Normalized coverage",
  featureAnnotation = "gene_name",
  just.group = "left",
  CDS = "darkred", ncRNA = "darkgreen",
  showFeatureId = TRUE,
  main = "Region 5", cex.main = 3
)

## Region 6
Gviz::plotTracks(
  list(ax, accDT3, accDT5, geneTrack),
  from = reg_starts[6], to = reg_ends[6],
  showId = TRUE, geneSymbol = TRUE,
  fontcolor.feature = "darkblue",
  name = "Normalized coverage",
  featureAnnotation = "gene_name",
  just.group = "left",
  CDS = "darkred", ncRNA = "darkgreen",
  showFeatureId = TRUE,
  main = "Region 6", cex.main = 3
)

## Region 7
Gviz::plotTracks(
  list(ax, accDT3, accDT5, geneTrack),
  from = reg_starts[7], to = reg_ends[7],
  showId = TRUE, geneSymbol = TRUE,
  fontcolor.feature = "darkblue",
  name = "Normalized coverage",
  featureAnnotation = "gene_name",
  just.group = "left",
  CDS = "darkred", ncRNA = "darkgreen",
  showFeatureId = TRUE,
  main = "Region 7", cex.main = 3
)

## Region 8
Gviz::plotTracks(
  list(ax, accDT3, accDT5, geneTrack),
  from = reg_starts[8], to = reg_ends[8],
  showId = TRUE, geneSymbol = TRUE,
  fontcolor.feature = "darkblue",
  name = "Normalized coverage",
  featureAnnotation = "gene_name",
  just.group = "left",
  CDS = "darkred", ncRNA = "darkgreen",
  showFeatureId = TRUE,
  main = "Region 8", cex.main = 3
)

## Region 9
Gviz::plotTracks(
  list(ax, accDT3, accDT5, geneTrack),
  from = reg_starts[9], to = reg_ends[9],
  showId = TRUE, geneSymbol = TRUE,
  fontcolor.feature = "darkblue",
  name = "Normalized coverage",
  featureAnnotation = "gene_name",
  just.group = "left",
  CDS = "darkred", ncRNA = "darkgreen",
  showFeatureId = TRUE,
  main = "Region 9", cex.main = 3
)

## Region 10
Gviz::plotTracks(
  list(ax, accDT3, accDT5, geneTrack),
  from = reg_starts[10], to = reg_ends[10],
  showId = TRUE, geneSymbol = TRUE,
  fontcolor.feature = "darkblue",
  name = "Normalized coverage",
  featureAnnotation = "gene_name",
  just.group = "left",
  CDS = "darkred", ncRNA = "darkgreen",
  showFeatureId = TRUE,
  main = "Region 10", cex.main = 3
)

## Highlight all regions in one canvas
ht <- Gviz::HighlightTrack(
  trackList = list(accDT3, accDT5, geneTrack),
  start     = reg_starts,
  end       = reg_ends
)
Gviz::plotTracks(
  list(ax, ht),
  showId = TRUE, geneSymbol = TRUE,
  just.group = "above",
  CDS = "darkred", ncRNA = "darkgreen"
)

