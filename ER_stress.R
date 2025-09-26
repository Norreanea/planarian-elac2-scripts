if (!require("pacman")) install.packages("pacman") 
#install required libraries (if needed) and load them
pacman::p_load(fdrtool,tximport,tximportData,edgeR,readr,tidyverse,RColorBrewer,pheatmap,DEGreport,ggplot2,ggrepel,"vsn",goeveg,zinbwave,genefilter,DESeq2,
               pheatmap,tidyverse,hrbrthemes,viridis,data.table,reshape,gplots,ggpubr,VennDiagram,limma,WGCNA,sm,network,venn,xlsx,openxlsx,glmGamPoi,Glimma ,
               dplyr,topGO,ALL,AnnotationDbi,org.Hs.eg.db,pathview,gage,gageData,plyr,strex,UpSetR,GOfuncR,gridExtra,grid,data.table,ggh4x,PCAtools,apeglm,sva,IHW)
library(ggplot2)
library(ggrepel)
#library(ggvenn)
pacman::p_load(ggvenn)
packageVersion("DESeq2")

library(DESeq2)
library(limma)
library(sva)
library(ggplot2)
library("IHW")
load("E:/Stringtie_anno/SM_anno/final/replicate4_5/local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23.RData")
#local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23

#Chech S-adenosylmethionine synthase SMESG000027587.1 SMESG000053202.1 SMESG000013953.1
load("E:/Stringtie_anno/SM_anno/final/replicate4_5/local_Wald_DESeq_corrected_dpa5_reduced_Elac15_WT52_GFP24.RData")
#local_Wald_DESeq_corrected_dpa5_reduced_Elac15_WT52_GFP24



local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23$all_genes


ER_stress_dpa3 <- local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23$all_genes
ER_stress_dpa3 <- ER_stress_dpa3[ER_stress_dpa3$gene %in% c("SMESG000044436.1",
                                                            "SMESG000058914.1",
                                                            "SMESG000058915.1",
                                                            "SMESG000054557.1",
                                                            "SMESG000046770.1",
                                                            "SMESG000021519.1",
                                                            "SMESG000064203.1",
                                                            "SMESG000068215.1",
                                                            "SMESG000010656.1",
                                                            "SMESG000054253.1",
                                                            "SMESG000017711.1",
                                                            "SMESG000061196.1",
                                                            "SMESG000048111.1",
                                                            "SMESG000065729.1"),]
ER_stress_dpa3 <- ER_stress_dpa3[,-3]
ER_stress_dpa3 <- unique(ER_stress_dpa3)

gene_map_ER <- data.frame(
  HS_gene = c("BLOC1S1", "SCARA3", "SCARA3", "PDGFRB", "Col6A1", "GALNT10", "XBP1", "BIP", "PDIA4", "HYOU1", "SEL1L", "EDEM1", "XBP1", "NQO1"),
  gene = c("SMESG000044436.1", "SMESG000058914.1", "SMESG000058915.1", "SMESG000054557.1", "SMESG000046770.1", "SMESG000021519.1", "SMESG000064203.1", "SMESG000068215.1", "SMESG000010656.1", "SMESG000054253.1", "SMESG000017711.1", "SMESG000061196.1", "SMESG000048111.1", "SMESG000065729.1"),
  stringsAsFactors = FALSE
)
ER_stress_dpa3 <- merge(ER_stress_dpa3, gene_map_ER, by = "gene", all.x = TRUE)


ER_stress_dpa5 <- local_Wald_DESeq_corrected_dpa5_reduced_Elac15_WT52_GFP24$all_genes
ER_stress_dpa5 <- ER_stress_dpa5[ER_stress_dpa5$gene %in% c("SMESG000044436.1",
                                                            "SMESG000058914.1",
                                                            "SMESG000058915.1",
                                                            "SMESG000054557.1",
                                                            "SMESG000046770.1",
                                                            "SMESG000021519.1",
                                                            "SMESG000064203.1",
                                                            "SMESG000068215.1",
                                                            "SMESG000010656.1",
                                                            "SMESG000054253.1",
                                                            "SMESG000017711.1",
                                                            "SMESG000061196.1",
                                                            "SMESG000048111.1",
                                                            "SMESG000065729.1"),]
ER_stress_dpa5 <- ER_stress_dpa5[,-3]
ER_stress_dpa5 <- unique(ER_stress_dpa5)
ER_stress_dpa5 <- merge(ER_stress_dpa5, gene_map_ER, by = "gene", all.x = TRUE)


ER_stress_dpa3$condition <- ifelse(ER_stress_dpa3$condition=="Elac_dpa3","ELAC2_KD",
                                    ifelse(ER_stress_dpa3$condition=="WT_dpa3","WT",
                                           ifelse(ER_stress_dpa3$condition=="GFP_dpa3","GFP_Mock","")))

ER_stress_dpa5$condition <- ifelse(ER_stress_dpa5$condition=="Elac_dpa5","ELAC2_KD",
                                   ifelse(ER_stress_dpa5$condition=="WT_dpa5","WT",
                                          ifelse(ER_stress_dpa5$condition=="GFP_dpa5","GFP_Mock","")))

ER_stress_dpa3$dpa <- "dpa3"
ER_stress_dpa5$dpa <- "dpa5"

ER_stress <- rbind(ER_stress_dpa3, ER_stress_dpa5)


library(ggplot2)
library(pheatmap)
library(dplyr)
# Prepare a matrix: rows=HS_gene, columns=sample; values=norm_counts
ER_stress <- as.data.frame(ER_stress)
mat <- ER_stress %>%
  dplyr::group_by(HS_gene, sample) %>%
  dplyr::summarise(count=mean(norm_counts), .groups="drop") %>%
  tidyr::pivot_wider(names_from=sample, values_from=count, values_fill=0) %>%
  tibble::column_to_rownames("HS_gene") %>%
  as.matrix()
# Create a vector mapping each sample to its condition
sample_conditions <- ER_stress %>%
  dplyr::select(sample, condition) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("sample")

# Make sure the rownames of annotation_col match colnames(mat)
annotation_col <- sample_conditions[colnames(mat), , drop = FALSE]

pheatmap(mat, annotation_col = annotation_col, 
         main = "Expression of ER stress/UPR genes (not DE in ELAC2 KD)")


ER_stress_summary <- ER_stress %>%
  dplyr::group_by(HS_gene, condition, dpa) %>%
  dplyr::summarise(mean_exp = mean(norm_counts), .groups="drop")


ER_stress_summary$sample <- paste(ER_stress_summary$condition, ER_stress_summary$dpa, sep = "_")

ER_stress_summary_long <- ER_stress_summary[,c(1,4,5)] %>%
  dplyr::select(HS_gene, mean_exp, sample) %>%
  tidyr::pivot_wider(names_from = sample, values_from = mean_exp, values_fill = 0) %>% 
  tibble::column_to_rownames("HS_gene") %>%
  as.matrix()

#?pheatmap
#scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
pheatmap(ER_stress_summary_long,
         main = "Expression of ER stress/UPR genes (not DE in ELAC2 KD)")



library(reshape2)   # used only if you prefer the ggplot version

#reorder the columns 
desired_order <- c("ELAC2_KD_dpa3", "ELAC2_KD_dpa5",
                   "GFP_Mock_dpa3", "GFP_Mock_dpa5",
                   "WT_dpa3",       "WT_dpa5")

mat_ord <- ER_stress_summary_long[ , desired_order]

#heat-map with pheatmap (no dendrograms) 
pheatmap(
  mat_ord,
  cluster_rows = FALSE,          # no row dendrogram
  cluster_cols = FALSE,          # no column dendrogram
  color = colorRampPalette(c("grey99", "red"))(100),   # white→red gradient
  #main  = "Expression of ER-stress / UPR genes (not DE in ELAC2 KD)",
  border_color = NA            
)
#Expression of ER-stress / UPR genes (not DE in ELAC2 KD)

ggplot(ER_stress_summary, aes(x = condition, y = HS_gene, size = mean_exp,color = dpa)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(title = "ER Stress/UPR Genes – No DE in ELAC2 KD", 
       x = "Condition", y = "Gene (HS)")



ggplot(ER_stress, aes(x = dpa, y = norm_counts, fill = condition)) +
  geom_boxplot() +
  facet_wrap(~HS_gene, scales = "free_y") +
  theme_bw() +
  labs(title = "ER stress/UPR genes - No DE in ELAC2 KD ",
       y = "Normalized counts")

ER_stress %>%
  dplyr::group_by(HS_gene, condition) %>%
  dplyr::summarise(mean_exp = mean(norm_counts)) %>%
  ggplot(aes(x = mean_exp, y = HS_gene, color = condition)) +
  geom_point(size = 3) +
  geom_segment(aes(xend = 0, yend = HS_gene), linetype="dotted", color="gray") +
  theme_bw() +
  labs(title = "ER stress/UPR genes expression (all not DE)", x = "Mean normalized expression")

ggplot(ER_stress, aes(x = sample, y = HS_gene, size = norm_counts, fill = condition)) +
  geom_point(shape = 21, color = "black") +
  theme_bw() +
  labs(title = "Coverage and expression of ER stress/UPR genes", y = "HS gene")

library(gt)
#install.packages("gt")
head(ER_stress)
ER_stress_noSM <- data.frame(gene=rep(NA, 8),
                             sample=rep(NA, 8),
                             norm_counts=rep(NA, 8),
                             condition=rep(NA, 8),
                             HS_gene=c("PMP22","MANF","TXNIP","HERPUD1","CHOP",
                                       "GADD34",
                                       "NOXA",
                                       "HMOX1"),
                             
                             dpa=rep(NA, 8)
                             )
ER_stress_all_ER <- rbind(ER_stress, ER_stress_noSM)
coverage_tbl <- ER_stress_all_ER %>%
  dplyr::group_by(HS_gene) %>%
  dplyr::summarise(
    "SM gene detected" = ifelse(any(!is.na(gene)), "✔", "✘"),
    "DE in ELAC2 KD" = "✘"
  ) %>% gt()
coverage_tbl


ggplot(ER_stress, aes(x = condition, y = norm_counts, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  facet_wrap(~HS_gene, scales = "free_y", ncol = 4) +
  theme_classic() +
  labs(title = "", x="",y="Normalized counts (median of ratios)")


library(UpSetR)
gene_lists <- list(
  Reference = unique(ER_stress_all_ER$HS_gene),
  Measured = unique(ER_stress$HS_gene)
)
upset(fromList(gene_lists))



ER_stress_summary
ER_stress_summary$sample <- paste(ER_stress_summary$condition, ER_stress_summary$dpa, sep = "_")
mat_mean <- ER_stress_summary %>%
  dplyr::group_by(HS_gene, sample) %>%
  dplyr::summarise(count=mean_exp, .groups="drop") %>%
  tidyr::pivot_wider(names_from=sample, values_from=count, values_fill=0) %>%
  tibble::column_to_rownames("HS_gene") %>%
  as.matrix()

ER_stress_summary$sample <- paste(ER_stress_summary$condition, ER_stress_summary$dpa, sep = "_")

mat_mean <- mat_mean[, order(annotation_col$condition, annotation_col$dpa)]
annotation_col <- annotation_col[order(annotation_col$condition, annotation_col$dpa), ]
pheatmap(
  mat_mean,
  annotation_col = annotation_col,
  cluster_cols = FALSE,  # disables dendrogram for columns (X)
  main = "Expression of ER stress/UPR genes (not DE in ELAC2 KD)"
)

unique(annotation_col$dpa)
annotation_col$dpa[is.na(annotation_col$dpa)] <- "NA"
dpa_levels <- unique(annotation_col$dpa)
# Get colors for unique dpa (at least 3 needed for Set2)
dpa_colors <- c( "#E78AC3" ,"#FFD92F")
names(dpa_colors) <- dpa_levels

ann_colors = list(
  condition = c("ELAC2_KD"="firebrick", "GFP_Mock"="darkgreen", "WT"="dodgerblue"),
  dpa = dpa_colors
)


pheatmap(
  mat_mean,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_cols = FALSE,
  main = "Expression of ER stress/UPR genes (not DE in ELAC2 KD)"
)
