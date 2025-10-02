############################################################
# qPCR KD validation + eye phenotype (%)
############################################################

# ---------- Packages ----------
# install.packages(c("readxl","dplyr","ggplot2","ggpubr","rstatix","car"))

library(readxl)    
library(dplyr)     
library(ggplot2)   
library(ggpubr)    
library(rstatix)   
library(car)       

# ---------- Input locations ----------
base_raw <- "https://raw.githubusercontent.com/Norreanea/planarian-elac2-scripts/main/data"
xlsx_url <- file.path(base_raw, "data.xlsx")

# Download to a temp file (readxl cannot read directly from a URL):
xlsx_tmp <- tempfile(fileext = ".xlsx")
download.file(xlsx_url, destfile = xlsx_tmp, mode = "wb")

# - Sheet 1: ELAC2_KD (qPCR fold-change)
# - Sheet 2: ELAC2_obs (counts of abnormalities per Gene (condition))
xlsx_tmp="D:/Elac2/final_results/tables/data.xlsx"
ELAC2_KD  <- readxl::read_excel(xlsx_tmp, sheet = 1)
ELAC2_obs <- readxl::read_excel(xlsx_tmp, sheet = 2)

# ---------- ELAC2_KD (qPCR fold-change) ----------
# Make sure numeric and factors are correct
ELAC2_KD <- ELAC2_KD %>%
  dplyr::mutate(
    FC       = as.numeric(FC),
    dpa      = as.factor(dpa),
    Condition = factor(Condition, levels = c("ELAC2","GFP","WT"),
                       labels = c("ELAC2 KD", "GFP mock", "WT"))
  )

# ---------- ASSUMPTION TESTS + STATS (per dpa) ----------
#   1) Fit one-way ANOVA model FC ~ Condition
#   2) Shapiro–Wilk on model residuals (normality)
#   3) Levene’s test (homogeneity)
#   4) If both OK -> ANOVA + Tukey HSD
#      Else       -> Kruskal–Wallis + Dunn (BH adjusted)

# Helper to run per-dpa stats with fallback
run_stats_per_dpa <- function(df) {
  # Fit model
  fit <- lm(FC ~ Condition, data = df)
  # Residual normality
  shap <- shapiro_test(residuals(fit)) %>%
    dplyr::mutate(test = "Shapiro-Wilk (residuals)")
  # Homogeneity
  lev  <- levene_test(FC ~ Condition, data = df) %>%
    dplyr::mutate(test = "Levene")
  # Decide parametric vs non-parametric
  normal_ok <- shap$p.value > 0.05
  hom_ok    <- lev$p > 0.05
  if (normal_ok && hom_ok) {
    # Parametric path: ANOVA + Tukey
    anova_res <- anova_test(data = df, dv = FC, wid = NULL, between = Condition)
    tukey     <- df %>% tukey_hsd(FC ~ Condition) %>% adjust_pvalue() %>% add_significance("p.adj")
    list(
      method     = "ANOVA + Tukey HSD",
      shapiro    = shap,
      levene     = lev,
      omnibus    = anova_res,
      posthoc    = tukey
    )
  } else {
    # Non-parametric path: Kruskal–Wallis + Dunn
    kw   <- kruskal_test(FC ~ Condition, data = df)
    dunn <- dunn_test(data = df, FC ~ Condition, p.adjust.method = "BH") %>%
      dplyr::mutate(p.adj.signif = cut(p.adj,
                                breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                labels = c("****","***","**","*","ns")))
    list(
      method     = "Kruskal–Wallis + Dunn (BH)",
      shapiro    = shap,
      levene     = lev,
      omnibus    = kw,
      posthoc    = dunn
    )
  }
}

# Run per dpa and bind post-hoc for plotting labels
stats_by_dpa <- ELAC2_KD %>%
  dplyr::group_split(dpa) %>%
  setNames(levels(ELAC2_KD$dpa)) %>%
  lapply(run_stats_per_dpa)

# Collect post-hoc results with dpa tag
collect_posthoc <- function(dpa_level, stats_list) {
  ph <- stats_list$posthoc
  ph$dpa <- dpa_level
  ph
}
posthoc_all <- bind_rows(mapply(collect_posthoc,
                                names(stats_by_dpa),
                                stats_by_dpa,
                                SIMPLIFY = FALSE))


if (!all(c("group1","group2") %in% names(posthoc_all))) {
  # just a safeguard
  posthoc_all <- posthoc_all %>%
    rename(group1 = group1, group2 = group2)
}

# y-positions for annotation 
y_positions <- c(1.2, 1.4, 1.6, 1.8) 
posthoc_all <- posthoc_all %>%
  dplyr::group_by(dpa) %>%
  dplyr::mutate(y.position = y_positions[seq_len(n())]) %>%
  ungroup()

# ---------- qPCR boxplot with per-dpa post-hoc ----------
(p_fig2a <- ggboxplot(ELAC2_KD, x = "Condition", y = "FC",
                     width = 0.25, size = 0.8,
                     facet.by = "dpa", short.panel.labs = TRUE,
                     bxp.errorbar = TRUE) +
  stat_pvalue_manual(posthoc_all,
                     label = "p.adj.signif",
                     hide.ns = TRUE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0),
             pch = 21, aes(fill = Condition), size = 2.8, show.legend = TRUE) +
  theme_bw(base_size = 20) +
  theme(
    text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    legend.position = "none"
  ) +
  labs(y = "Fold change") +
  scale_x_discrete(labels = c("ELAC2 KD","GFP mock","WT")) +
  scale_y_continuous(limits = c(0, 1.8)))

# print(p_fig2a)
################################################################################
# ---------------------- ELAC2 phenotype --------------
# Calculate %
N_PER_CELL <- 10  # number of worms

# Harmonize Gene labels for the legend 
ELAC2_obs <- ELAC2_obs %>%
  dplyr::mutate(
    Gene = case_when(
      Gene %in% c("ELAC","ELAC2") ~ "ELAC2 KD",
      Gene == "GFP"               ~ "GFP mock",
      Gene == "WT"                ~ "WT",
      TRUE                        ~ Gene
    )
  )

# Select only tail part
ELAC2_obs <- ELAC2_obs[ELAC2_obs$Body_part=="tail",]
# percent = 100 * sum(Photoreceptor_abnormalities) / (N_PER_CELL * number_of_rows_in_group)
ELAC2_obs_pct <- ELAC2_obs %>%
  dplyr::group_by(Gene, dpa) %>%
  dplyr::summarise(
    n_rows = n(),
    abnormal_sum = sum(Photoreceptor_abnormalities, na.rm = TRUE),
    total_examined = N_PER_CELL * n_rows,
    pct_abnormal = 100 * abnormal_sum / total_examined,
    .groups = "drop"
  ) %>%
  # ensure dpa numeric for continuous x-axis
  dplyr::mutate(dpa = as.numeric(dpa))

all_groups <- c("WT","GFP mock","ELAC2 KD")
ELAC2_obs_pct <- ELAC2_obs_pct %>%
  # add dpa=0 rows (0%) for all groups
  bind_rows(tibble(Gene = all_groups, dpa = 1, pct_abnormal = 0))%>%
  bind_rows(tibble(Gene = all_groups, dpa = 2, pct_abnormal = 0))%>%
  bind_rows(tibble(Gene = all_groups, dpa = 3, pct_abnormal = 0))%>%
  bind_rows(tibble(Gene = all_groups, dpa = 4, pct_abnormal = 0))

as.data.frame(ELAC2_obs_pct)
unique(ELAC2_obs_pct$Gene)
# ---------- eye abnormalities over time (%) ----------
(p_fig2b <- ggplot(ELAC2_obs_pct, aes(x = dpa, y = pct_abnormal, color = Gene, group = Gene)) +
   #geom_line(na.rm = TRUE) +
   #geom_point(na.rm = TRUE) +
   # optional smoother (comment out if you prefer raw lines only)
   stat_smooth(se = FALSE, method = "loess", span = 0.3, alpha = 1, na.rm = TRUE,size=4) +
   theme_classic(base_size = 22) +
   labs(y = "Worms with eye \nabnormalities (%)", x = "dpa") +
   scale_y_continuous(limits = c(0, 100)) +
   scale_x_continuous(breaks = 0:dpa_max, limits = c(0, dpa_max)) +
   theme(
     axis.text.x       = element_text(size = 18),
     legend.title      = element_blank(),
     legend.text       = element_text(size = 18),
     panel.background  = element_blank(),
     panel.grid.major  = element_blank(),
     legend.position   = c(0.18, 0.78),
     legend.justification = "left"
   ))

# print(p_fig2b)

