############################################################
# Supplementary Figure S5 — ncRNA length profiles (S5A)
# and a schematic of the dominant 5′ tRNA-half (S5B)
# RData source for S5A is not provided
############################################################

## install if missing
if (!requireNamespace("readxl",  quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr",   quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr",   quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("scales",  quietly = TRUE)) install.packages("scales")
if (!requireNamespace("RRNA",    quietly = TRUE)) install.packages("RRNA")


## -------------------------------------------------------
## S5A: ncRNA general length distribution 
## -------------------------------------------------------
# Read and normalize column 

load("G:/Elac2/smallRNA/sncRNAlist_filtered_ELAC_filtered_sets_table.RData") #Huge object
tbl_len$length=nchar(tbl_len$seq)
#tbl_len <- readxl::read_excel(xlsx_path, sheet = sheet_lengths)

# Harmonize 
if (!"seq" %in% names(tbl_len) && "Sequence" %in% names(tbl_len)) {
  tbl_len <- dplyr::rename(tbl_len, seq = "Sequence")
}
if (!"correct_RNA_type" %in% names(tbl_len) && "RNA_type" %in% names(tbl_len)) {
  tbl_len <- dplyr::rename(tbl_len, correct_RNA_type = "RNA_type")
}

# Compute length
tbl_len <- tbl_len |>
  dplyr::mutate(length = nchar(.data$seq)) |>
  dplyr::filter(!.data$correct_RNA_type %in% c("multiple_hits", "genome")) |>
  dplyr::mutate(correct_RNA_type = dplyr::recode(.data$correct_RNA_type,
                                                 "other" = "other fragments",
                                                 .default = .data$correct_RNA_type))

# Factor order 
lvl_types <- c("piRNA","mRNA fragments","lncRNA fragments",
               "tRNA fragments","snoRNA fragments","rRNA fragments",
               "snRNA fragments","miRNA","other fragments")
tbl_len$correct_RNA_type <- factor(tbl_len$correct_RNA_type, levels = lvl_types, ordered = TRUE)

# Colors
pal_types <- c(
  "piRNA"="#1b9e77","mRNA fragments"="lightblue","lncRNA fragments"="#7570b3",
  "tRNA fragments"="salmon","snoRNA fragments"="#66a61e","rRNA fragments"="grey",
  "snRNA fragments"="#a6761d","miRNA"="#e6ab02","other fragments"="#1f78b4"
)

# Plot S5A
S5A <- ggplot2::ggplot(tbl_len, ggplot2::aes(x = length, color = correct_RNA_type, fill = correct_RNA_type)) +
  ggplot2::geom_histogram(alpha = 0.6, binwidth = 1, position = "identity") +
  ggplot2::scale_fill_manual(values = pal_types, drop = FALSE) +
  ggplot2::scale_color_manual(values = pal_types, drop = FALSE) +
  ggplot2::scale_y_continuous(labels = scales::label_scientific()) +
  ggplot2::scale_x_continuous(breaks = seq(15, 115, by = 20)) +
  ggplot2::facet_wrap(~ correct_RNA_type, scales = "free") +
  ggplot2::labs(x = "Length", y = "Number of RNA species", subtitle = "ncRNA length distributions") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    legend.position = "none",
    panel.spacing   = grid::unit(0.1, "lines"),
    strip.text.x    = ggplot2::element_text(size = 12, hjust = 0.5, family = "Times"),
    axis.title.x    = ggplot2::element_text(size = 12, family = "Times"),
    axis.title.y    = ggplot2::element_text(size = 12, family = "Times"),
    text            = ggplot2::element_text(family = "Times")
  )

print(S5A)

## -------------------------------------------------------
## S5B: 5′ tRNA-half secondary structure panel (RRNA::RNAPlot)
## -------------------------------------------------------
dom <- data.frame(
  Structure="(((((((..((((.......)))).(((((.......)))))....(((((.......)))))))))))).",
  rnacentral_Seq="GCAUCGGUGGUUCAGUGGUAGAAUACUCGCCUGCCACGCGGGCGGCCCGGGUUAGAUUCCCGGCCGAUGCA",
  Sequence="GCAUCGGUGGUUCAGUGGUAGAAUACUCGCCU",
  start=1,
  end=32,
  tRF_ID="5'-tRNA half Gly-GCC"
  
)
#dom <- dom[1, , drop = FALSE]
# Convert dot-bracket if structure contains '<' '>' 
struct_raw <- as.character(dom$Structure)
struct_db  <- gsub(">", ")", gsub("<", "(", struct_raw), fixed = TRUE)

seq_full   <- as.character(dom$rnacentral_Seq)
seq_range  <- c(as.integer(dom$start), as.integer(dom$end))
label_txt  <- paste("tRF", as.character(dom$Sequence))

# Build CT and coordinates
ct_obj   <- RRNA::makeCt(struct_db, seq_full)
coords   <- RRNA::ct2coord(ct_obj)
rangesDF <- data.frame(
  min  = seq_range[1],
  max  = seq_range[2],
  col  = 2,
  desc = label_txt
)

# Rotate the layout a bit for readability
rotate_coords <- function(coords, angle_deg = 30) {
  rad <- angle_deg * pi / 180
  R   <- matrix(c(cos(rad), -sin(rad), sin(rad), cos(rad)), ncol = 2)
  out <- coords
  xy  <- as.matrix(out[, c("x","y")])
  out[, c("x","y")] <- t(R %*% t(xy))
  out
}
coords_rot <- rotate_coords(coords, angle_deg = 45)

# Plot S5B
RRNA::RNAPlot(coords_rot,
              nt      = TRUE,
              ranges  = rangesDF,
              labTF   = TRUE,
              main    = as.character(dom$tRF_ID),
              pointSize = 3,
              tsize     = 1)
