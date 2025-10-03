############################################################
# Figure 7 — Smed ELAC2 KD remodels the sncRNA landscape
# Panels:
# A  RNA class composition 
# B  Differential sncRNA species
# C  miRNAs changed by ELAC2 KD (3 dpa & 5 dpa; N/D labels)
# D  qPCR validation of selected miRNAs (3 dpa) 
# E  Differential tRNA fragments: by fragment class and amino acid
# F  qPCR of 5′ tiRNA Gly-GCC (3 dpa) 
############################################################

## Core packages 
if (!requireNamespace("readxl", quietly = TRUE))  install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE))   install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE))   install.packages("tidyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggpubr", quietly = TRUE))  install.packages("ggpubr")
if (!requireNamespace("rstatix", quietly = TRUE)) install.packages("rstatix")
if (!requireNamespace("ggh4x", quietly = TRUE))   install.packages("ggh4x")
if (!requireNamespace("scales", quietly = TRUE))  install.packages("scales")
if (!requireNamespace("plotrix", quietly = TRUE)) install.packages("plotrix")
library(ggplot2)


use_local <- TRUE 

if (use_local) {
  xlsx_path <- "D:/Elac2/final_results/tables/data.xlsx"
} else {
  gh_base   <- "https://raw.githubusercontent.com/Norreanea/planarian-elac2-scripts/main/data"
  xlsx_url  <- file.path(gh_base, "data.xlsx")
  xlsx_tmp  <- tempfile(fileext = ".xlsx")
  utils::download.file(xlsx_url, destfile = xlsx_tmp, mode = "wb", quiet = TRUE)
  xlsx_path <- xlsx_tmp
}

## ------------------------------------------------------------------
## SHEET MAP (adjust names if your workbook uses different labels)
## ------------------------------------------------------------------
sheet_sncRNA_DGE      <- "snRNA_DEG_no_rRNA"    
sheet_all_RNA_diff    <- "all_snRNA"           
sheet_ncRNA_validation<- "ncRNA_validation"      
sheet_tRF_classified  <- "tRF_classified"         
sheet_RNA_distribution<- "snRNA_distribution"       

## ------------------------------------------------------------------
## Load data from data.xlsx by sheet name
## ------------------------------------------------------------------
not_rRNA <- readxl::read_excel(xlsx_path, sheet = sheet_sncRNA_DGE)
tRNA_old <- readxl::read_excel(xlsx_path, sheet = sheet_all_RNA_diff)
ncRNA_val <- readxl::read_excel(xlsx_path, sheet = sheet_ncRNA_validation)
tRF_class <- readxl::read_excel(xlsx_path, sheet = sheet_tRF_classified)
RNA_dist  <- readxl::read_excel(xlsx_path, sheet = sheet_RNA_distribution)

## ------------------------------------------------------------------
## Panel pre-processing
## ------------------------------------------------------------------


if (!"dpa" %in% names(not_rRNA) && "set" %in% names(not_rRNA)) {
  not_rRNA$dpa <- ifelse(not_rRNA$set %in% c("Elac_vs_WT_dpa3","GFP_vs_WT_dpa3"), "dpa3", "dpa5")
}

# (E) tRNA fragment helpers 
if ("snRNA_type" %in% names(not_rRNA)) {
  extract_tRNA_and_amino <- function(x) {
    parts <- strsplit(x, "-", fixed = TRUE)[[1]]
    if (length(parts) < 3) return(c(NA, NA))
    c(paste(parts[1:2], collapse = "-"), paste(parts[3:length(parts)], collapse = "-"))
  }
  idx <- which(not_rRNA$RNA_type == "tRNA fragments" & !is.na(not_rRNA$snRNA_type))
  if (length(idx)) {
    tmp <- t(vapply(not_rRNA$snRNA_type[idx], extract_tRNA_and_amino, FUN.VALUE = character(2)))
    not_rRNA$tRNA_structure <- NA_character_
    not_rRNA$amino          <- NA_character_
    not_rRNA$tRNA_structure[idx] <- tmp[,1]
    not_rRNA$amino[idx]          <- tmp[,2]
  }
}

## ------------------------------------------------------------------
## (F) qPCR: 5′ tiRNA Gly-GCC (ELAC2 KD vs WT), Tukey on 2 groups → t-test
## ------------------------------------------------------------------
tRNA_val <- ncRNA_val |>
  dplyr::filter(!startsWith(.data$ncRNA, "Sme") & .data$ncRNA %in% c("5` tRNA half Gly-GCC-1","5` tRNA half Gly-GCC")) |>
  dplyr::mutate(
    FC        = as.numeric(.data$FC),
    Condition = dplyr::recode(.data$Condition, "ELAC2" = "ELAC2 KD"),
    ncRNA     = "5′ tRNA half Gly-GCC"
  )

# two-sided t-test
stat.test_tRNA <- rstatix::t_test(tRNA_val, FC ~ Condition, detailed = TRUE) |>
  rstatix::add_significance("p") |>
  dplyr::mutate(y.position = max(tRNA_val$FC, na.rm = TRUE) * 1.1)

F_fig7 <-
  ggpubr::ggboxplot(tRNA_val, x = "Condition", y = "FC",
                    width = 0.25, size = 0.2, fill = "Condition") +
  ggpubr::stat_pvalue_manual(stat.test_tRNA, label = "p.signif",
                             y.position = "y.position", tip.length = 0.01) +
  geom_point(position = position_jitterdodge(jitter.width = 2, dodge.width = 0),
             pch = 21, aes(fill = Condition), size = 2, show.legend = FALSE) +
  labs(y = "Fold change", x = NULL, subtitle = "5′ tiRNA Gly-GCC (qPCR, 3 dpa)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        panel.grid = element_blank())

## ------------------------------------------------------------------
## (D) qPCR: selected miRNAs (3 dpa), two-sided t-test per miRNA
## ------------------------------------------------------------------
miRNA_val <- ncRNA_val |>
  dplyr::filter(startsWith(.data$ncRNA, "Sme")) |>
  dplyr::mutate(
    FC        = as.numeric(.data$FC),
    Condition = dplyr::recode(.data$Condition, "ELAC2" = "ELAC2 KD"),
    ncRNA     = gsub("Sme-miR-", "sme-miR-", .data$ncRNA)
  )

stat.test_miRNA <- miRNA_val |>
  dplyr::group_by(ncRNA) |>
  rstatix::t_test(FC ~ Condition, detailed = TRUE) |>
  rstatix::add_significance("p") |>
  dplyr::ungroup() |>
  dplyr::group_by(ncRNA) |>
  dplyr::mutate(y.position = max(miRNA_val$FC[miRNA_val$ncRNA == unique(ncRNA)], na.rm = TRUE) * 1.1) |>
  dplyr::ungroup()

D_fig7 <-
  ggpubr::ggboxplot(miRNA_val, x = "Condition", y = "FC",
                    width = 0.25, size = 0.2, fill = "Condition",
                    facet.by = "ncRNA", ncol = 1, short.panel.labs = TRUE) +
  ggpubr::stat_pvalue_manual(stat.test_miRNA, label = "p.signif",
                             y.position = "y.position", tip.length = 0.01) +
  geom_point(position = position_jitterdodge(jitter.width = 2, dodge.width = 0),
             pch = 21, aes(fill = Condition), size = 2, show.legend = FALSE) +
  labs(y = "Fold change", x = NULL, subtitle = "miRNA qPCR (3 dpa)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        panel.grid = element_blank())

## ------------------------------------------------------------------
## (E) Differential tRNA fragments by class and amino acid
## ------------------------------------------------------------------
only_tRNA <- tRF_class |>
  dplyr::transmute(
    tRNA_structure = dplyr::recode(.data$class,
                                   "whole-mature"    = "mature tRNA",
                                   "whole-pre-tRNA"  = "pre-tRNA",
                                   "5′-tRNA-half"    = "5′ tiRNA",
                                   "3′-tRNA-half"    = "3′ tiRNA",
                                   "5′-tRF"          = "5′ tRF",
                                   "3′-tRF"          = "3′ tRF",
                                   .default          = .data$class
    ),
    amino        = .data$final_family,
    localizarion = dplyr::recode(.data$localizarion, "mitochondria" = "mtDNA", .default = .data$localizarion),
    dpa          = .data$dpa,
    log2FoldChange = .data$log2FoldChange
  ) |>
  dplyr::filter(!.data$tRNA_structure %in% c("mature tRNA", "pre-tRNA")) |>
  dplyr::mutate(
    # nicer capitalization for amino 
    amino = sub("([A-Z])([A-Z]{2})", "\\1\\L\\2", .data$amino, perl = TRUE)
  )

only_tRNA$tRNA_structure <- ifelse(only_tRNA$amino == "ALA-TGC", "i-tRF", only_tRNA$tRNA_structure)

E_fig7 <-
  ggplot2::ggplot(only_tRNA, ggplot2::aes(x = tRNA_structure, y = amino)) +
  ggplot2::geom_jitter(ggplot2::aes(color = log2FoldChange, shape = dpa),
                       size = 3, width = 0.2, height = 0.05) +
  ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  ggh4x::facet_grid2(localizarion ~ tRNA_structure, scales = "free", space = "free_y") +
  ggplot2::labs(x = NULL, y = NULL, color = "log2FC", shape = "dpa",
                subtitle = "Differential tRNA fragments by fragment class and amino acid") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    strip.background   = element_rect(fill = "white"),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom= element_blank(),
    panel.spacing      = unit(0, "lines"),
    strip.text.x       = element_text(size = 10)
  )

## ------------------------------------------------------------------
## (C) miRNAs changed in ELAC2 KD; add N/D labels, points by dpa & log2FC
## ------------------------------------------------------------------
# "sme-miR-190a-3p"  - D     
# "sme-mir-2153_precursor" - N
# "sme-miR-8b-3p"          
# "sme-miR-87c-5p"         
# "sme-miR-36a-3p"  -N       
# "sme-miR-125b-5p"   -D     
# "sme-miR-9a-5p"   - D       
# "sme-miR-71a-5p"  - N      
# "sme-miR-96b-5p"   -D      
# "sme-miR-96a-5p"  - D
miRNA <- not_rRNA |>
  dplyr::filter(.data$RNA_type == "miRNA") |>
  dplyr::mutate(
    ND = dplyr::case_when(
      .data$snRNA_type %in% c("sme-mir-2153_precursor","sme-miR-36a-3p","sme-miR-71a-5p") ~ "N",
      .data$snRNA_type %in% c("sme-miR-8b-3p","sme-miR-87c-5p")                            ~ "",
      TRUE ~ "D"
    ),
    snRNA_type_new = paste(.data$snRNA_type, .data$ND)
  )

C_fig7 <-
  ggplot2::ggplot(miRNA, ggplot2::aes(x = snRNA_type_new, y = "")) +
  ggplot2::geom_jitter(ggplot2::aes(color = log2FoldChange, shape = dpa),
                       size = 3, width = 0.4, height = 0.1) +
  ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  ggh4x::facet_grid2(~ snRNA_type_new, scales = "free", space = "free_y") +
  ggplot2::labs(x = NULL, y = NULL, color = "log2FC", shape = "dpa",
                subtitle = "miRNAs differentially accumulated in ELAC2 KD") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    strip.background   = element_rect(fill = "white"),
    axis.text.x.bottom = element_text(angle = 90, hjust = 1),
    axis.ticks.x.bottom= element_blank(),
    panel.spacing      = unit(0, "lines"),
    strip.text.x       = element_text(size = 10)
  )

## ------------------------------------------------------------------
## (A) RNA class composition
## ------------------------------------------------------------------
RNA_dist$`RNA type` <- factor(
  RNA_dist$`RNA type`,
  levels = c("piRNA","mRNA fragments","lncRNA fragments","tRNA fragments",
             "snoRNA fragments","rRNA fragments","snRNA fragments","miRNA","other fragments"),
  ordered = TRUE
)

# label helper for facet gene naming
RNA_dist$gene <- dplyr::recode(RNA_dist$gene,
                               "ELAC2" = "ELAC2\nKD",
                               "GFP"   = "GFP\nmock",
                               .default = RNA_dist$gene)

custom_cols_A <- c(
  "piRNA"="#1b9e77","mRNA fragments"="lightblue","lncRNA fragments"="#7570b3",
  "tRNA fragments"="salmon","snoRNA fragments"="#66a61e","rRNA fragments"="grey",
  "snRNA fragments"="#a6761d","miRNA"="#e6ab02","other fragments"="#1f78b4"
)

A_fig7 <-
  ggplot2::ggplot(RNA_dist |> dplyr::arrange(gene, dplyr::desc(`RNA type`)) |>
                    dplyr::group_by(gene) |>
                    dplyr::mutate(text_y = cumsum(perc) - perc/2) |>
                    dplyr::ungroup(),
                  ggplot2::aes(x = "", y = perc, fill = `RNA type`)) +
  ggplot2::geom_bar(stat = "identity", width = 1.5) +
  ggplot2::scale_fill_manual(values = custom_cols_A) +
  ggplot2::geom_label(ggplot2::aes(label = paste0(perc, "%"), y = text_y),
                      show.legend = FALSE, label.size = 0.1) +
  ggh4x::facet_nested(. ~ dpa + gene) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.title = element_text(size = 12),
    axis.ticks.x = element_blank()
  ) +
  ggplot2::labs(y = "Total normalized read count (%)", x = NULL, fill = "RNA species",
                subtitle = "RNA class composition")

## ------------------------------------------------------------------
## (B) Pyramid: number of sncRNA species with accumulation vs loss 
## ------------------------------------------------------------------
pyr_df <- not_rRNA |>
  dplyr::group_by(RNA_type, expression, dpa) |>
  dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
  tidyr::pivot_wider(names_from = expression, values_from = count, values_fill = 0) |>
  dplyr::rename(accumulation = "accumulation", loss = "loss") |>
  dplyr::mutate(loss = -loss)

rna_order <- c("lncRNA fragments","miRNA","other fragments","snoRNA fragments","tRNA fragments")

pyr_dpa <- function(df, which_dpa) {
  out <- df |>
    dplyr::filter(dpa == which_dpa) |>
    dplyr::select(RNA_type, accumulation, loss) |>
    dplyr::arrange(factor(RNA_type, levels = rna_order)) |>
    dplyr::mutate(loss = abs(loss)) |>
    tidyr::complete(RNA_type = rna_order, fill = list(accumulation = 0, loss = 0))
  out
}

pyramid_dpa3 <- pyr_dpa(pyr_df, "dpa3")
pyramid_dpa5 <- pyr_dpa(pyr_df, "dpa5")

## Complex base-graphics pyramid panel (
graphics::grid.newpage()
graphics::layout(matrix(c(1, 2, 3), nrow = 1), widths = c(3, 1.2, 3))

# Left (dpa3)
graphics::par(family = "Times", mar = c(0, 0, 0, 0), ps = 12)
max_val_loss3 <- max(pyramid_dpa3$loss)
max_val_acc3  <- max(pyramid_dpa3$accumulation)
plotrix::pyramid.plot(
  rx      = pyramid_dpa3$accumulation,
  lx      = pyramid_dpa3$loss,
  labels  = rep("", nrow(pyramid_dpa3)),
  lxcol   = c(scales::alpha("#7570b3", 0.6), scales::alpha("#e6ab02", 0.6),
              scales::alpha("#1f78b4", 0.6), scales::alpha("#66a61e", 0.6),
              scales::alpha("salmon", 0.6)),
  rxcol   = c(scales::alpha("#7570b3", 0.6), scales::alpha("#e6ab02", 0.6),
              scales::alpha("#1f78b4", 0.6), scales::alpha("#66a61e", 0.6),
              scales::alpha("salmon", 0.6)),
  main    = "dpa3",
  top.labels = c("Loss","","Accumulation"),
  unit    = "Number of RNA species",
  labelcex= 0.8,
  show.values = TRUE,
  xlim    = c(-max_val_loss3, max_val_acc3),
  laxlab  = seq(0, max_val_loss3, by = 10),
  raxlab  = seq(0, max_val_acc3,  by = 10),
  ndig    = 0
)

# Middle (RNA type labels)
graphics::plot.new()
lbls <- c("lncRNA \nfragments","miRNA","other \nfragments","snoRNA \nfragments","tRNA \nfragments")
graphics::text(0.5, seq(0.1, 0.9, length.out = length(lbls)), labels = lbls, cex = 1)

# Right (dpa5)
graphics::par(mar = c(0, 0, 0, 0))
max_val_loss5 <- max(pyramid_dpa5$loss) + 3
max_val_acc5  <- max(pyramid_dpa5$accumulation) + 3
plotrix::pyramid.plot(
  rx      = pyramid_dpa5$accumulation,
  lx      = pyramid_dpa5$loss,
  labels  = rep("", nrow(pyramid_dpa5)),
  lxcol   = c(scales::alpha("#7570b3", 0.6), scales::alpha("#e6ab02", 0.6),
              scales::alpha("#1f78b4", 0.6), scales::alpha("#66a61e", 0.6),
              scales::alpha("salmon", 0.6)),
  rxcol   = c(scales::alpha("#7570b3", 0.6), scales::alpha("#e6ab02", 0.6),
              scales::alpha("#1f78b4", 0.6), scales::alpha("#66a61e", 0.6),
              scales::alpha("salmon", 0.6)),
  main    = "dpa5",
  top.labels = c("Loss","","Accumulation"),
  gap     = 0.4,
  unit    = "Number of RNA species",
  labelcex= 0.8,
  show.values = TRUE,
  xlim    = c(-max_val_loss5, max_val_acc5),
  laxlab  = seq(0, max_val_loss5, by = 5),
  raxlab  = seq(0, max_val_acc5,  by = 5),
  ndig    = 0
)

## ------------------------------------------------------------------
## Arrange panels (C, D, E, F)
## ------------------------------------------------------------------
ggpubr::ggarrange(C_fig7, D_fig7, E_fig7, F_fig7,
                  ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")

