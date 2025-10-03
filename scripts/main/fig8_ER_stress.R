############################################################
# Figure 8 – ER-stress
############################################################
# install.packages(c("readr","dplyr","tidyr","tibble","pheatmap","ggplot2"))
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(ggplot2)

## Read data from
gh_base   <- "https://raw.githubusercontent.com/Norreanea/planarian-elac2-scripts/main/data"
dpa3  <- file.path(gh_base, "norm_counts_3dpa.xlsx")  # columns: gene, condition, sample, norm_counts
dpa5  <- file.path(gh_base, "norm_counts_3dpa.xlsx")


## ER-stress gene list (Smed gene IDs and human labels for plotting)
er_ids <- c("SMESG000044436.1","SMESG000058914.1","SMESG000058915.1","SMESG000054557.1",
            "SMESG000046770.1","SMESG000021519.1","SMESG000064203.1","SMESG000068215.1",
            "SMESG000010656.1","SMESG000054253.1","SMESG000017711.1","SMESG000061196.1",
            "SMESG000048111.1","SMESG000065729.1")

gene_map_ER <- data.frame(
  gene   = er_ids,
  HS_gene = c("BLOC1S1","SCARA3","SCARA3","PDGFRB","Col6A1","GALNT10","XBP1","BIP",
              "PDIA4","HYOU1","SEL1L","EDEM1","XBP1","NQO1"),
  stringsAsFactors = FALSE
)

## Load dpa3/dpa5 data
read_er_xlsx <- function(url_xlsx, er_ids) {
  tf <- base::tempfile(fileext = ".xlsx")
  utils::download.file(url_xlsx, tf, mode = "wb", quiet = TRUE)
  df <- readxl::read_xlsx(tf)
  
  # Expect columns: gene, condition, sample, norm_counts
  df |>
    dplyr::filter(.data$gene %in% er_ids) |>
    dplyr::select(dplyr::any_of(c("gene","condition","sample","norm_counts"))) |>
    dplyr::mutate(
      gene        = as.character(.data$gene),
      condition   = as.character(.data$condition),
      sample      = as.character(.data$sample),
      norm_counts = as.numeric(.data$norm_counts)
    ) |>
    dplyr::distinct()
}

er_dpa3 <- read_er_xlsx(dpa3, er_ids = er_ids) |>
  dplyr::mutate(dpa = "dpa3")

er_dpa5 <- read_er_xlsx(dpa5, er_ids = er_ids) |>
  dplyr::mutate(dpa = "dpa5")

## Clean condition labels and tag dpa 
normalize_cond <- function(x) {
  dplyr::case_when(
    x %in% c("Elac_dpa3","Elac_dpa5") ~ "ELAC2_KD",
    x %in% c("WT_dpa3","WT_dpa5")     ~ "WT",
    x %in% c("GFP_dpa3","GFP_dpa5")   ~ "GFP_Mock",
    TRUE                              ~ x
  )
}

er_dpa3 <- er_dpa3 |>
  dplyr::mutate(condition = normalize_cond(.data$condition),
                dpa       = "dpa3")
er_dpa5 <- er_dpa5 |>
  dplyr::mutate(condition = normalize_cond(.data$condition),
                dpa       = "dpa5")

## Merge, map human gene symbols, and keep only necessary columns
ER_stress <- dplyr::bind_rows(er_dpa3, er_dpa5) |>
  dplyr::inner_join(gene_map_ER, by = "gene") |>
  dplyr::mutate(
    norm_counts = as.numeric(.data$norm_counts),
    sample      = as.character(.data$sample),
    condition   = as.character(.data$condition),
    dpa         = as.character(.data$dpa),
    HS_gene     = as.character(.data$HS_gene)
  )

## Heatmap 2 – means per (condition × dpa) 
ER_stress_summary <- ER_stress |>
  dplyr::group_by(.data$HS_gene, .data$condition, .data$dpa) |>
  dplyr::summarise(mean_exp = mean(.data$norm_counts, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(sample = paste(.data$condition, .data$dpa, sep = "_"))

mat_cond_dpa <- ER_stress_summary |>
  dplyr::select(.data$HS_gene, .data$mean_exp, .data$sample) |>
  tidyr::pivot_wider(names_from = .data$sample, values_from = .data$mean_exp, values_fill = 0) |>
  tibble::column_to_rownames("HS_gene") |>
  as.matrix()

# Optional: order columns in a biologically sensible sequence
desired_order <- c("ELAC2_KD_dpa3","ELAC2_KD_dpa5",
                   "GFP_Mock_dpa3","GFP_Mock_dpa5",
                   "WT_dpa3","WT_dpa5")
# Keep only columns that exist; reorder them
desired_order <- intersect(desired_order, colnames(mat_cond_dpa))
mat_ord <- mat_cond_dpa[, desired_order, drop = FALSE]

# Per-condition×dpa means heatmap (no clustering, white→red)
pheatmap::pheatmap(
  mat_ord,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = grDevices::colorRampPalette(c("grey99", "red"))(100),
  main  = "Expression of ER stress genes (not DE in ELAC2 KD)",
  border_color = NA
)
