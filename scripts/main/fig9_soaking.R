############################################################
# Figure 9 – 5′ tiRNA Gly-GCC soaking experiments
############################################################

## Packages 
# install.packages(c("readxl","dplyr","ggplot2","ggpubr","rstatix","car","emmeans","multcomp","stringr","tidyr"))
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)    
library(rstatix)   
library(car)       
library(emmeans)   
library(multcomp)  
library(stringr)
library(tidyr)

# ---------- Input locations ----------
base_raw <- "https://raw.githubusercontent.com/Norreanea/planarian-elac2-scripts/main/data"
xlsx_url <- file.path(base_raw, "data.xlsx")

# Download to a temp file (readxl cannot read directly from a URL):
xlsx_tmp <- tempfile(fileext = ".xlsx")
utils::download.file(xlsx_url, destfile = xlsx_tmp, mode = "wb")

# For local testing keep your path; comment out if using GitHub raw file above
xlsx_tmp <- "D:/Elac2/final_results/tables/data.xlsx"

tRF   <- readxl::read_excel(xlsx_tmp, sheet = 3)
genes <- readxl::read_excel(xlsx_tmp, sheet = 4)

# ensure numeric
tRF$Expression   <- as.numeric(tRF$Expression)
genes$Expression <- as.numeric(genes$Expression)

## Harmonize condition names 
harmonize_condition <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "5'?tiRNAGlyGCC", "5'tiRNA Gly-GCC")
  x <- stringr::str_replace_all(x, "5` tiRNA half[[:space:]]*GLY-GCC", "5'tiRNA Gly-GCC")
  x <- stringr::str_replace_all(x, "5` tiRNA[[:space:]]*GLY-GCC", "5'tiRNA Gly-GCC")
  x <- stringr::str_replace_all(x, "5'\\s*tRNA.*Gly.*GCC", "5'tiRNA Gly-GCC")
  x <- stringr::str_replace_all(x, "^ELAC2\\s*KD$", "ELAC2 KD")
  x <- stringr::str_replace_all(x, "^GFP\\s*Mock$", "GFP Mock")
  x <- stringr::str_replace_all(x, "^WT$", "WT")
  x <- stringr::str_replace_all(x, "planarian\\s*water", "planarian water")
  x <- stringr::str_replace_all(x, "scrambled", "scrambled")
  x
}

tRF$Condition   <- harmonize_condition(tRF$Condition)
genes$Condition <- harmonize_condition(genes$Condition)

## Factor orders for x-axis 
cond_levels <- c(
  "ELAC2 KD 5'tiRNA Gly-GCC","ELAC2 KD planarian water","ELAC2 KD scrambled",
  "GFP Mock 5'tiRNA Gly-GCC","GFP Mock planarian water","GFP Mock scrambled",
  "WT 5'tiRNA Gly-GCC","WT planarian water","WT scrambled"
)
tRF$Condition   <- factor(tRF$Condition, levels = cond_levels)
genes$Condition <- factor(genes$Condition, levels = cond_levels)

## Colors 
pal <- c(
  "ELAC2 KD 5'tiRNA Gly-GCC"="#F8766D","ELAC2 KD planarian water"="#F8766D","ELAC2 KD scrambled"="#F8766D",
  "GFP Mock 5'tiRNA Gly-GCC"="#00BA38","GFP Mock planarian water"="#00BA38","GFP Mock scrambled"="#00BA38",
  "WT 5'tiRNA Gly-GCC"="#619CFF","WT planarian water"="#619CFF","WT scrambled"="#619CFF"
)

# -------------------------------------------------------------------
#  Assumption helpers
# -------------------------------------------------------------------
# assumption_checks <- function(df){
#   fit  <- stats::lm(Expression ~ Condition, data = df)
#   shap <- rstatix::shapiro_test(stats::residuals(fit))
#   lev  <- rstatix::levene_test(Expression ~ Condition, data = df)
#   list(shapiro_p = shap$p.value, levene_p = lev$p, normal = shap$p.value > 0.05, homo = lev$p > 0.05)
# }
# helper: map p-values to significance stars
.p_to_stars <- function(p) {
  dplyr::case_when(
    p <= 1e-4 ~ "****",
    p <= 1e-3 ~ "***",
    p <= 1e-2 ~ "**",
    p <= 5e-2 ~ "*",
    TRUE      ~ "ns"
  )
}
# lmer + Tukey 
do_lmm_glht_tukey <- function(df){
  if (!"Replicate" %in% names(df)) stop("Replicate column required for lmer.")
  suppressPackageStartupMessages(require(lme4))
  m  <- lme4::lmer(Expression ~ Condition + (1 | Replicate), data = df, REML = TRUE)
  
  gl <- multcomp::glht(m, linfct = multcomp::mcp(Condition = "Tukey"))
  sm <- summary(gl)  # (Adjusted p values reported -- single-step method)
  
  td <- broom::tidy(gl) |>
    tidyr::separate(contrast, into = c("group1","group2"), sep = " - ") |>
    dplyr::mutate(
      p.adj        = .data$adj.p.value,
      p.adj.signif = .p_to_stars(.data$adj.p.value)
    )
  
  list(model = m, posthoc = td, method = "lmer + Tukey (multcomp single-step)")
}


# do_kw_dunn <- function(df){
#   kw   <- rstatix::kruskal_test(df, Expression ~ Condition)
#   dunn <- rstatix::dunn_test(df, Expression ~ Condition, p.adjust.method = "BH") |>
#     dplyr::mutate(p.adj.signif = dplyr::case_when(
#       p.adj <= 1e-4 ~ "****",
#       p.adj <= 1e-3 ~ "***",
#       p.adj <= 1e-2 ~ "**",
#       p.adj <= 5e-2 ~ "*",
#       TRUE ~ "ns"
#     ))
#   list(omnibus = kw, posthoc = dunn)
# }
# 
# run_stats <- function(df){
#   A <- assumption_checks(df)
#   if (A$normal && A$homo) {
#     S <- do_lmm_tukey(df); S$method <- "LMM + Tukey (parametric)"
#   } else {
#     S <- do_kw_dunn(df);   S$method <- "Kruskal–Wallis + Dunn (BH)"
#   }
#   S$assumptions <- A
#   S
# }

# -------------------------------------------------------------------
# ANOVA + Tukey  for group4/group5
# Justification: ANOVA is reasonably robust to mild violations of normality and 
# homoscedasticity—especially with comparable group sizes—so it’s acceptable to 
# use ANOVA + Tukey for Group 4/5, while I`m keeping the mixed-effects model for the 
# gene-level analysis.
# -------------------------------------------------------------------
do_anova_tukey <- function(df){
  aov_fit <- stats::aov(Expression ~ Condition, data = df)
  emm     <- emmeans::emmeans(aov_fit, ~ Condition)
  # FIX: use contrast(..., "pairwise") instead of pairs()
  tuk     <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") |>
    summary(infer = TRUE)
  
  post    <- tuk |>
    dplyr::mutate(contrast = as.character(contrast)) |>
    tidyr::separate(contrast, into = c("group1","group2"), sep = " - ") |>
    dplyr::mutate(p.adj = p.value,
                  p.adj.signif = dplyr::case_when(
                    p.adj <= 1e-4 ~ "****",
                    p.adj <= 1e-3 ~ "***",
                    p.adj <= 1e-2 ~ "**",
                    p.adj <= 5e-2 ~ "*",
                    TRUE ~ "ns"
                  ))
  list(model = aov_fit, posthoc = post, method = "ANOVA + Tukey (forced)")
}

##  tRF-level plots (Groups 4 and 5) – now ANOVA + Tukey ----------------------
plot_tRF_group <- function(dat, group_name, colors = pal){
  df <- dat |>
    dplyr::filter(Group == group_name) |>
    tidyr::drop_na(Condition, Expression)
  if (nrow(df) < 3) stop("Not enough data in ", group_name)
  
  # Force ANOVA + Tukey
  S <- do_anova_tukey(df)
  post <- S$posthoc |> dplyr::filter(p.adj.signif != "ns")
  #post <- S$posthoc
  # bracket positions
  ymax <- max(df$Expression, na.rm = TRUE)
  if (nrow(post) > 0) {
    post <- post |>
      dplyr::mutate(y.position = ymax * (1.05 + 0.07 * (dplyr::row_number()-1)))
  }
  
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = Condition, y = Expression, fill = Condition)) +
    ggplot2::geom_boxplot(width = .6, outlier.shape = NA) +
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .15), size = 2, alpha = .7) +
    { if (nrow(post) > 0)
      ggpubr::stat_pvalue_manual(post,
                                 label      = "p.adj.signif",
                                 y.position = "y.position",
                                 tip.length = 0.01,
                                 xmin       = "group1", xmax = "group2",
                                 inherit.aes= FALSE)
      else NULL } +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, .05))) +
    ggplot2::labs(subtitle = paste0(group_name, "  (", S$method, ")"),
                  y = "Relative 5′-tiRNA Gly-GCC accumulation", x = NULL) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "none")
  gg
}

(p_group4 <- plot_tRF_group(tRF, "group4"))
(p_group5 <- plot_tRF_group(tRF, "group5"))

################################################################################
## Gene-level bar/SE plot for the PTP target (SMESG000048842.1)  ----

target_id <- "SMESG000048842.1"

gene_sum <- genes |>
  dplyr::filter(Target == target_id) |>
  dplyr::group_by(Condition, Target) |>
  dplyr::summarise(mean_exp = stats::median(Expression, na.rm = TRUE),
                   sd_exp   = stats::sd(Expression, na.rm = TRUE),
                   n        = dplyr::n(),
                   se_exp   = sd_exp / sqrt(n),
                   .groups  = "drop")

gene_df <- genes |> dplyr::filter(Target == target_id)

# IMPORTANT: match your old path
S_gene <- do_lmm_glht_tukey(gene_df)

# Keep ONLY these 3 comparisons (order doesn’t matter, both "A - B" and "B - A" possible)
ctrl <- "ELAC2 KD 5'tiRNA Gly-GCC"
keepers <- list(
  c(ctrl, "ELAC2 KD planarian water"),
  c(ctrl, "ELAC2 KD scrambled"),
  c(ctrl, "WT planarian water")
)

post_g <- S_gene$posthoc |>
  dplyr::rowwise() |>
  dplyr::filter(
    any(c(group1, group2) == ctrl) &&
      paste(sort(c(group1, group2)), collapse = " | ") %in%
      sapply(keepers, function(k) paste(sort(k), collapse = " | "))
  ) |>
  dplyr::ungroup()

# Make sure direction is ctrl (group1) vs other (group2) for bracket mapping
post_g <- post_g |>
  dplyr::mutate(
    g1 = ifelse(group1 == ctrl, group1, group2),
    g2 = ifelse(group1 == ctrl, group2, group1),
    group1 = g1, group2 = g2
  ) |>
  dplyr::select(-g1, -g2)

# y-positions for the 3 brackets
ymax_g <- max(gene_sum$mean_exp + gene_sum$se_exp, na.rm = TRUE)
post_g <- post_g |>
  dplyr::mutate(y.position = ymax_g * (1.05 + 0.07 * (dplyr::row_number()-1)))
# Color by background group (unchanged)
gene_sum <- gene_sum |>
  dplyr::mutate(GeneGrp = dplyr::case_when(
    stringr::str_starts(as.character(Condition), "ELAC2 KD") ~ "ELAC2 KD",
    stringr::str_starts(as.character(Condition), "GFP Mock") ~ "GFP Mock",
    stringr::str_starts(as.character(Condition), "WT")       ~ "WT",
    TRUE ~ "Other"
  ))

(p_gene <- ggplot2::ggplot(gene_sum, ggplot2::aes(x = Condition, y = mean_exp, fill = GeneGrp)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.8) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_exp - se_exp, ymax = mean_exp + se_exp),
                         width = 0.2, position = ggplot2::position_dodge(width = 0.8)) +
  ggpubr::stat_pvalue_manual(
    post_g,
    label      = "p.adj.signif",
    y.position = "y.position",
    tip.length = 0.01,
    xmin       = "group1", xmax = "group2",
    inherit.aes= FALSE
  ) +
  ggplot2::scale_fill_manual(values = c("ELAC2 KD"="#F8766D","GFP Mock"="#00BA38","WT"="#619CFF")) +
  ggplot2::labs(title = target_id,
                y = "Relative gene expression", x = NULL,
                subtitle = "") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))


# Print plots
p_group4
p_group5
p_gene
