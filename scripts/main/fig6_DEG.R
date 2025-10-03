#--------------------------load DESeq input objects from GitHub-----------------

base_raw <- "https://raw.githubusercontent.com/Norreanea/planarian-elac2-scripts/main/data" 


load_RData_from_github <- function(relpath) {
  url  <- paste0(base_raw, "/", relpath)
  tf   <- tempfile(fileext = ".RData")
  download.file(url, tf, mode = "wb", quiet = TRUE)
  load(tf, envir = .GlobalEnv)
}


load_RData_from_github("all_stringtie_dds_dpa3.RData")
load_RData_from_github("repo/all_stringtie_dds_dpa5.RData")

# sanity check 
stopifnot(exists("all_stringtie_dds_dpa3"), exists("all_stringtie_dds_dpa5"))


#------------------------------libraries----------------------------------------

#packageVersion("DESeq2") #1.44.20
#install required libraries (if needed) and load them
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse, dplyr, tidyr, data.table, ggplot2, ggrepel, hrbrthemes, viridis, ggpubr, ggtext,
  readr, openxlsx,
  rtracklayer, GenomicFeatures, Gviz, txdbmaker,
  DESeq2, IHW, sva, limma,
  topGO,
  ggh4x, cowplot, patchwork
)

#-----------------------------------main----------------------------------------
my_DESeq <- function(mycountData, subset_option = "full", test = "Wald", fitType = "local", 
                     betaPrior = FALSE, sfType = "ratio", cooksCutoff = TRUE, reduced = NULL, minReplicatesForReplace = 7,
                     independentFiltering = TRUE, alpha = 0.05, pAdjustMethod = "BH", 
                     parallel = FALSE, minmu = 0.5) {
  # Subset DESeqDataSet according to the option
  mycountData_sub=mycountData
  
  if (subset_option != "full") {
    cols_to_keep <- switch(subset_option,
                           "dpa3" = grep("3_batch1$|3_batch2$|3_batch3$|L13$|L23$|L33$|L43$|L53$", colnames(mycountData), value = TRUE),
                           "dpa5" = grep("5_batch1$|5_batch2$|5_batch3$|L15$|L25$|L35$|L45$|L55$", colnames(mycountData), value = TRUE),
                           "reduced" = grep("Elac_L23_batch1|Elac_L23_batch2|Elac_L15_batch1|Elac_L15_batch2|GFP_L33_batch1|GFP_L33_batch2|GFP_L25_batch1|GFP_L25_batch2", colnames(mycountData), value = TRUE, invert = TRUE),
                           "dpa3_reduced_Elac23_WT23_GFP23" = grep("3_batch1$|3_batch2$|3_batch3$|L13$|L23$|L33$|L43$|L53$", 
                                                                   grep("^Elac_L23|^Elac_L33|^WT_L23|^WT_L33|^GFP_L23|^GFP_L33", 
                                                                        colnames(mycountData), value = TRUE, invert = TRUE), value = TRUE),
                           "dpa5_reduced_Elac12_WT12_GFP24" = grep("5_batch1$|5_batch2$|5_batch3$|L15$|L25$|L35$|L45$|L55$", 
                                                                   grep("^Elac_L25|^Elac_L15|^WT_L25|^WT_L15|^GFP_L25|^GFP_L45", 
                                                                        colnames(mycountData), value = TRUE, invert = TRUE), value = TRUE),
                           "dpa3_reduced_Elac15_WT53_GFP23" = grep("3_batch1$|3_batch2$|3_batch3$|L13$|L23$|L33$|L43$|L53$", 
                                                                   grep("^Elac_L13|^Elac_L53|^WT_L53|^WT_L33|^GFP_L23|^GFP_L33", 
                                                                        colnames(mycountData), value = TRUE, invert = TRUE), value = TRUE),
                           "dpa5_reduced_Elac15_WT52_GFP24" = grep("5_batch1$|5_batch2$|5_batch3$|L15$|L25$|L35$|L45$|L55$", 
                                                                   grep("^Elac_L15|^Elac_L55|^WT_L55|^WT_L25|^GFP_L25|^GFP_L45", 
                                                                        colnames(mycountData), value = TRUE, invert = TRUE), value = TRUE),
                           "dpa3_reduced_Elac25_WT53_GFP23" = grep("3_batch1$|3_batch2$|3_batch3$|L13$|L23$|L33$|L43$|L53$", 
                                                                   grep("^Elac_L23|^Elac_L53|^WT_L53|^WT_L33|^GFP_L23|^GFP_L33", 
                                                                        colnames(mycountData), value = TRUE, invert = TRUE), value = TRUE),
                           "dpa5_reduced_Elac15_WT52_GFP23" = grep("5_batch1$|5_batch2$|5_batch3$|L15$|L25$|L35$|L45$|L55$", 
                                                                   grep("^Elac_L15|^Elac_L55|^WT_L55|^WT_L25|^GFP_L25|^GFP_L35", 
                                                                        colnames(mycountData), value = TRUE, invert = TRUE), value = TRUE),
                           "dpa3_reduced_Elac21_WT53_GFP23" = grep("3_batch1$|3_batch2$|3_batch3$|L13$|L23$|L33$|L43$|L53$", 
                                                                   grep("^Elac_L23|^Elac_L13|^WT_L53|^WT_L33|^GFP_L23|^GFP_L33", 
                                                                        colnames(mycountData), value = TRUE, invert = TRUE), value = TRUE),
                           
                           
    )
    #dpa3_reduced_Elac21_WT53_GFP23, dpa3_reduced_Elac25_WT53_GFP23
    #dpa5_reduced_Elac15_WT52_GFP23, dpa5_reduced_Elac15_WT52_GFP24 
    
    mycountData_sub <- mycountData[, cols_to_keep]
  }
  run=colnames(mycountData_sub)
  stringtie_samples=as.data.frame(run)
  if (all(grepl("batch", colnames(mycountData_sub)))) {
    (stringtie_samples$sample=substr(stringtie_samples$run,1,nchar(stringtie_samples$run)-7))
    (stringtie_samples$gene=substr(stringtie_samples$run,1,nchar(stringtie_samples$run)-11))
    (stringtie_samples$replicate=paste("rep",substr(stringtie_samples$sample,nchar(stringtie_samples$sample)-1,nchar(stringtie_samples$sample)-1),sep=""))
    (stringtie_samples$platform=ifelse(stringtie_samples$replicate %in% c("rep4","rep5"),"NovaSeq","NextSeq"))
    (stringtie_samples$dpa=paste("dpa",substr(stringtie_samples$sample,nchar(stringtie_samples$sample),nchar(stringtie_samples$sample)),sep=""))
    (stringtie_samples$batch=substr(stringtie_samples$run,nchar(stringtie_samples$run),nchar(stringtie_samples$run)))
    (stringtie_samples$condition=as.factor(paste(stringtie_samples$gene,stringtie_samples$dpa,sep="_")))
  } else {
    (stringtie_samples$sample=stringtie_samples$run)
    (stringtie_samples$gene=substr(stringtie_samples$run,1,nchar(stringtie_samples$run)-4))
    (stringtie_samples$replicate=paste("rep",substr(stringtie_samples$sample,nchar(stringtie_samples$sample)-1,nchar(stringtie_samples$sample)-1),sep=""))
    (stringtie_samples$platform=ifelse(stringtie_samples$replicate %in% c("rep4","rep5"),"NovaSeq","NextSeq"))
    (stringtie_samples$dpa=paste("dpa",substr(stringtie_samples$sample,nchar(stringtie_samples$sample),nchar(stringtie_samples$sample)),sep=""))
    (stringtie_samples$condition=as.factor(paste(stringtie_samples$gene,stringtie_samples$dpa,sep="_")))
    
  }
  
  # all_stringtie_dds=DESeqDataSetFromMatrix(countData = mycountData_sub,
  #                                          colData = stringtie_samples,
  #                                          design = ~ condition + platform)
  all_stringtie_dds=DESeqDataSetFromMatrix(countData = mycountData_sub,
                                           colData = stringtie_samples,
                                           design = ~ platform + condition)
  # all_stringtie_dds=DESeqDataSetFromMatrix(countData = mycountData_sub,
  #                                          colData = stringtie_samples,
  #                                          design = ~ condition)
  ddsColl_stringtie = collapseReplicates(all_stringtie_dds, all_stringtie_dds$sample, all_stringtie_dds$run)
  
  
  
  # Obtain a list of conditions dynamically
  conditions <- unique(ddsColl_stringtie$condition)
  dds_sub <- estimateSizeFactors(ddsColl_stringtie,type="ratio") 
  # Initialize a vector to keep track of genes to keep
  keep_non_zero_genes <- rep(FALSE, nrow(dds_sub))
  keep_high_expr_genes <- rep(FALSE, nrow(dds_sub))
  
  # Loop through each condition to check for counts criteria
  for (cond in conditions) {
    # Get samples for this condition
    samples_in_cond <- colnames(dds_sub)[dds_sub$condition == cond]
    # Apply the filtering criteria
    keep_for_this_group <- rowSums(counts(dds_sub, normalized = TRUE)[, samples_in_cond] >= 10) >= 2
    keep_high_expr_genes <- keep_high_expr_genes | keep_for_this_group
    keep_for_non_zero <- rowSums(counts(dds_sub, normalized = TRUE)[, samples_in_cond] >= 1) >= 2
    keep_non_zero_genes <- keep_non_zero_genes | keep_for_non_zero
  }
  
  # Filter the DESeqDataSet
  dds_filtered <- dds_sub[keep_high_expr_genes, ]
  dds_non_zero <- dds_sub[keep_non_zero_genes, ]
  
  dds_filtered <- estimateSizeFactors(dds_filtered,type="ratio") 
  
  # Basic statistics after subsetting
  print(paste("Number of non zero genes in ", subset_option,": ",nrow(dds_non_zero)))
  print(paste("Number genes with at least 10 counts in ", subset_option,": ",nrow(dds_filtered)))
  
  
  # Prepare for DESeq analysis
  if (test == "LRT" && !is.null(reduced)) {
    dds_new <- DESeq(dds_sub, test = test, fitType = fitType, betaPrior = betaPrior,
                     full = design(dds_sub), reduced = reduced, quiet = FALSE,
                     minReplicatesForReplace = minReplicatesForReplace, sfType = sfType,
                     parallel = parallel)
  } else {
    dds_new <- DESeq(dds_sub, test = test, fitType = fitType, betaPrior = betaPrior,
                     full = design(dds_sub), quiet = FALSE, minReplicatesForReplace = minReplicatesForReplace,
                     sfType = sfType, parallel = parallel)
  }
  

  contrasts <- list(
    GFP_vs_WT_dpa3 = c("condition", "GFP_dpa3", "WT_dpa3"),
    Elac_vs_GFP_dpa3 = c("condition", "Elac_dpa3", "GFP_dpa3"),
    Elac_vs_WT_dpa3 = c("condition", "Elac_dpa3", "WT_dpa3"),
    GFP_vs_WT_dpa5 = c("condition", "GFP_dpa5", "WT_dpa5"),
    Elac_vs_GFP_dpa5 = c("condition", "Elac_dpa5", "GFP_dpa5"),
    Elac_vs_WT_dpa5 = c("condition", "Elac_dpa5", "WT_dpa5"),
    WT_dpa3_vs_dpa5 = c("condition", "WT_dpa3", "WT_dpa5"),
    GFP_dpa3_vs_dpa5 = c("condition", "GFP_dpa3", "GFP_dpa5"),
    Elac_dpa3_vs_dpa5 = c("condition", "Elac_dpa3", "Elac_dpa5")
  )
  
  contrasts_filtered <- Filter(function(x) {
    # Check if both the second and third elements of each contrast are in conditions
    all(c(x[2], x[3]) %in% conditions)
  }, contrasts)
  
  # Results storage
  results_list <- list()
  # Collect all p-values from all contrasts
  all_pvalues <- c()
  
  for (contrast_name in names(contrasts_filtered)) {
    res <- results(dds_new, contrast = contrasts_filtered[[contrast_name]], cooksCutoff = cooksCutoff,
                   independentFiltering = independentFiltering, alpha = alpha,
                   pAdjustMethod = pAdjustMethod, parallel = parallel, minmu = minmu)
    # Filter results
    res <- res[!is.na(res$padj) & !is.na(res$pvalue), ]
    res$set <- contrast_name
    res$gene <- rownames(res)
    results_list[[contrast_name]] <- res
    all_pvalues <- c(all_pvalues, res$pvalue)
  }
  
  # Adjust p-values across all contrasts
  adjusted_pvalues <- p.adjust(all_pvalues, method = "BH")
  
  # Reassign adjusted p-values back to the results
  index <- 1
  for (contrast_name in names(results_list)) {
    res <- results_list[[contrast_name]]
    res$padj <- adjusted_pvalues[index:(index + nrow(res) - 1)]
    results_list[[contrast_name]] <- res
    index <- index + nrow(res)
  }
  # 
  
  # Combine and return results
  combined_results <- do.call(rbind, results_list)
  # Filtering for padj <= 0.05
  pvalue_good <- subset(combined_results, padj <= 0.05)
  
  #Add normalized counts
  gene_norm_count=as.data.frame(counts(dds_new, normalized=TRUE))
  gene_norm_count$gene=rownames(gene_norm_count)
  #DGE_set_correct=pvalue_good
  pvalue_good=pvalue_good[pvalue_good$gene!="SMESG000037458.1",]
  DGE_raw_set_with_norm_counts=merge(as.data.frame(pvalue_good),gene_norm_count,by="gene",all.x=TRUE)
  combined_results_with_norm_counts=merge(as.data.frame(combined_results),gene_norm_count,by="gene",all.x=TRUE)
  #colnames(DGE_raw_set_with_norm_counts)
  
  #Check if results are reliable:
  test_DGE=DGE_raw_set_with_norm_counts
  ELAC3=which((startsWith(colnames(test_DGE),"Elac")) & (endsWith(colnames(test_DGE),"3")))
  ELAC5=which((startsWith(colnames(test_DGE),"Elac")) & (endsWith(colnames(test_DGE),"5")))
  WT3=which((startsWith(colnames(test_DGE),"WT")) & (endsWith(colnames(test_DGE),"3")))
  WT5=which((startsWith(colnames(test_DGE),"WT")) & (endsWith(colnames(test_DGE),"5")))
  GFP3=which((startsWith(colnames(test_DGE),"GFP")) & (endsWith(colnames(test_DGE),"3")))
  GFP5=which((startsWith(colnames(test_DGE),"GFP")) & (endsWith(colnames(test_DGE),"5")))
  
  # Initialize new columns for condition1 and condition2 indices
  test_DGE$columns_condition1 <- vector("list", nrow(test_DGE))
  test_DGE$columns_condition2 <- vector("list", nrow(test_DGE))
  
  # # Loop through each row
  for (i in 1:nrow(test_DGE)) {
    # Use switch to determine the columns for condition1 and condition2
    test_DGE$columns_condition1[[i]] <- switch(test_DGE$set[i],
                                               "GFP_vs_WT_dpa3" = GFP3,
                                               "Elac_vs_GFP_dpa3" = ELAC3,
                                               "GFP_dpa3_vs_dpa5" = GFP3,
                                               "GFP_vs_WT_dpa5" = GFP5,
                                               "Elac_vs_WT_dpa3" = ELAC3,
                                               "Elac_vs_GFP_dpa5" = ELAC5,
                                               "WT_dpa3_vs_dpa5" = WT3,
                                               "Elac_dpa3_vs_dpa5" = ELAC5,
                                               "Elac_vs_WT_dpa5" = ELAC5,
                                               NULL
    )
    
    test_DGE$columns_condition2[[i]] <- switch(test_DGE$set[i],
                                               "GFP_vs_WT_dpa3" = WT3,
                                               "Elac_vs_GFP_dpa3" = GFP3,
                                               "GFP_dpa3_vs_dpa5" = GFP5,
                                               "GFP_vs_WT_dpa5" = WT5,
                                               "Elac_vs_WT_dpa3" = WT3,
                                               "Elac_vs_GFP_dpa5" = GFP5,
                                               "WT_dpa3_vs_dpa5" = WT5,
                                               "Elac_dpa3_vs_dpa5" = ELAC5,
                                               "Elac_vs_WT_dpa5" = WT5,
                                               NULL
    )
    test_DGE$consistency[i]=ifelse((length(which(test_DGE[i,c(unlist(test_DGE$columns_condition1[[i]]))] >0))>=2|
                                      length(which(test_DGE[i,c(unlist(test_DGE$columns_condition2[[i]]))] >0))>=2),
                                   "reliable","not_reliable")
  }
  DGE_raw_set_with_norm_counts_consistensy_tested=test_DGE
  DGE_raw_set_with_norm_counts_consistensy_tested[DGE_raw_set_with_norm_counts_consistensy_tested$consistency=="not_reliable",]
  test_DGE=test_DGE[test_DGE$consistency=="reliable",]
  
  test_DGE=test_DGE[abs(test_DGE$log2FoldChange)>0.58,]
  test_DGE_GFP_WT_dpa3=test_DGE[test_DGE$set=="GFP_vs_WT_dpa3",]
  test_DGE_GFP_WT_dpa5=test_DGE[test_DGE$set=="GFP_vs_WT_dpa5",]
  #Include ELac_vs_GFP
  final_test_DGE_no_threshold=rbind(test_DGE[test_DGE$set=="Elac_vs_WT_dpa3" & test_DGE$log2FoldChange>=0 & (!(test_DGE$gene %in% 
                                                                                                                 test_DGE$gene[test_DGE$set=="GFP_vs_WT_dpa3"& 
                                                                                                                                 test_DGE$log2FoldChange>=0])|
                                                                                                               test_DGE$gene %in% 
                                                                                                               test_DGE$gene[test_DGE$set=="Elac_vs_GFP_dpa3"& 
                                                                                                                               test_DGE$log2FoldChange>=0]),] ,
                                    
                                    test_DGE[test_DGE$set=="Elac_vs_WT_dpa3" & test_DGE$log2FoldChange<0 & (!(test_DGE$gene %in% 
                                                                                                                test_DGE$gene[test_DGE$set=="GFP_vs_WT_dpa3"& 
                                                                                                                                test_DGE$log2FoldChange<0])|
                                                                                                              test_DGE$gene %in% 
                                                                                                              test_DGE$gene[test_DGE$set=="Elac_vs_GFP_dpa3"& 
                                                                                                                              test_DGE$log2FoldChange<0]),],
                                    test_DGE[test_DGE$set=="Elac_vs_WT_dpa5" & test_DGE$log2FoldChange>=0 & (!(test_DGE$gene %in% 
                                                                                                                 test_DGE$gene[test_DGE$set=="GFP_vs_WT_dpa5"& 
                                                                                                                                 test_DGE$log2FoldChange>=0])|
                                                                                                               test_DGE$gene %in% 
                                                                                                               test_DGE$gene[test_DGE$set=="Elac_vs_GFP_dpa5"& 
                                                                                                                               test_DGE$log2FoldChange>=0]),],
                                    test_DGE[test_DGE$set=="Elac_vs_WT_dpa5" & test_DGE$log2FoldChange<0 & (!(test_DGE$gene %in% 
                                                                                                                test_DGE$gene[test_DGE$set=="GFP_vs_WT_dpa5"& 
                                                                                                                                test_DGE$log2FoldChange<0])|
                                                                                                              test_DGE$gene %in% 
                                                                                                              test_DGE$gene[test_DGE$set=="Elac_vs_GFP_dpa5"& 
                                                                                                                              test_DGE$log2FoldChange<0]),],
                                    test_DGE[test_DGE$set=="Elac_vs_GFP_dpa3",],
                                    test_DGE[test_DGE$set=="Elac_vs_GFP_dpa5",]
  )
  
  #VENN diagram for filtered DEGs
  DGE_set_correct=test_DGE
  DGE_set_correct$gene_id=DGE_set_correct$gene
  library(ggvenn)
  #?ggvenn
  up_dpa3=list(
    Elac_vs_WT_dpa3=DGE_set_correct$gene_id[DGE_set_correct$log2FoldChange>=0 & DGE_set_correct$set=="Elac_vs_WT_dpa3"],
    Elac_vs_GFP_dpa3=DGE_set_correct$gene_id[DGE_set_correct$log2FoldChange>=0 & DGE_set_correct$set=="Elac_vs_GFP_dpa3"],
    GFP_vs_WT_dpa3=DGE_set_correct$gene_id[DGE_set_correct$log2FoldChange>=0 & DGE_set_correct$set=="GFP_vs_WT_dpa3"]
  )
  up_dpa3_venn_filtered=ggvenn(
    up_dpa3, 
    fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 4
  )
  up_dpa5=list(
    Elac_vs_WT_dpa5=DGE_set_correct$gene_id[DGE_set_correct$log2FoldChange>=0 & DGE_set_correct$set=="Elac_vs_WT_dpa5"],
    Elac_vs_GFP_dpa5=DGE_set_correct$gene_id[DGE_set_correct$log2FoldChange>=0 & DGE_set_correct$set=="Elac_vs_GFP_dpa5"],
    GFP_vs_WT_dpa5=DGE_set_correct$gene_id[DGE_set_correct$log2FoldChange>=0 & DGE_set_correct$set=="GFP_vs_WT_dpa5"]
  )
  up_dpa5_venn_filtered=ggvenn(
    up_dpa5, 
    fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 4
  )
  down_dpa3=list(
    Elac_vs_WT_dpa3=DGE_set_correct$gene_id[DGE_set_correct$log2FoldChange<0 & DGE_set_correct$set=="Elac_vs_WT_dpa3"],
    Elac_vs_GFP_dpa3=DGE_set_correct$gene_id[DGE_set_correct$log2FoldChange<0 & DGE_set_correct$set=="Elac_vs_GFP_dpa3"],
    GFP_vs_WT_dpa3=DGE_set_correct$gene_id[DGE_set_correct$log2FoldChange<0 & DGE_set_correct$set=="GFP_vs_WT_dpa3"]
  )
  down_dpa3_venn_filtered=ggvenn(
    down_dpa3, 
    fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 4
  )
  down_dpa5=list(
    Elac_vs_WT_dpa5=DGE_set_correct$gene_id[DGE_set_correct$log2FoldChange<0 & DGE_set_correct$set=="Elac_vs_WT_dpa5"],
    Elac_vs_GFP_dpa5=DGE_set_correct$gene_id[DGE_set_correct$log2FoldChange<0 & DGE_set_correct$set=="Elac_vs_GFP_dpa5"],
    GFP_vs_WT_dpa5=DGE_set_correct$gene_id[DGE_set_correct$log2FoldChange<0 & DGE_set_correct$set=="GFP_vs_WT_dpa5"]
  )
  down_dpa5_venn_filtered=ggvenn(
    down_dpa5, 
    fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 4
  )
  
  
  
  vsd <- varianceStabilizingTransformation(dds_sub)
  # Print statements
  cat("Number of all unique DEG with padj < 0.05:", length(unique(pvalue_good$gene)), "\n")
  # reliable_genes=length(unique(pvalue_good$gene[abs(pvalue_good$log2FoldChange)>1,]))
  # reliable_genes <- unique(test_DGE$gene[test_DGE$consistency == "reliable"])
  cat("Number of all unique DEG with padj < 0.05 &  abs(lfc) >1:", length(unique(pvalue_good$gene[abs(pvalue_good$log2FoldChange)>1])), "\n")
  cat("Number of all unique DEG with padj < 0.05 in GFP vs WT dpa3:", length(unique(test_DGE_GFP_WT_dpa3$gene)), "\n")
  cat("Number of all unique DEG with padj < 0.05 in GFP vs WT dpa5:", length(unique(test_DGE_GFP_WT_dpa5$gene)), "\n")
  cat("Number of common in both days for GFP vs WT:", sum(unique(test_DGE_GFP_WT_dpa3$gene) %in% unique(test_DGE_GFP_WT_dpa5$gene)), "\n")
  cat("Common in both days for GFP vs WT:", test_DGE_GFP_WT_dpa3$gene[unique(test_DGE_GFP_WT_dpa3$gene) %in% unique(test_DGE_GFP_WT_dpa5$gene)], "\n")
  cat("Number of unique DEG with padj < 0.05  & not present in GFP vs WT:", length(unique(final_test_DGE_no_threshold$gene)), "\n")
  # cat("Number of unique DEG with padj < 0.05 & abs(logFC) > 1 & not present in GFP vs WT:", 
  #     length(unique(final_test_DGE_no_threshold$gene[abs(final_test_DGE_no_threshold$log2FoldChange)>=1])), "\n")
  cat("Number of unique DEG with padj < 0.05 & abs(logFC) > 0.58 & not present in GFP vs WT:", 
      length(unique(final_test_DGE_no_threshold$gene[abs(final_test_DGE_no_threshold$log2FoldChange)>=0.58])), "\n")
  
  #Get normalized counts for final DEG
  sigOE=final_test_DGE_no_threshold[,c("gene","set",grep("_L", colnames(final_test_DGE_no_threshold), value = TRUE))]
  sigOE_long=gather(sigOE,sample,norm_counts,colnames(sigOE)[-c(1,2)],factor_key=TRUE)
  sigOE_long=merge(sigOE_long,colData(dds_new)[,c("sample","condition")])
  
  #Get normalized counts for all genes
  combined_results_with_counts=merge(as.data.frame(combined_results),gene_norm_count,by="gene",all.x=TRUE)
  combined_results_with_counts=combined_results_with_counts[,c("gene","set",grep("_L", colnames(combined_results_with_counts), value = TRUE))]
  combined_results_with_counts_long=gather(combined_results_with_counts,sample,norm_counts,colnames(combined_results_with_counts)[-c(1,2)],factor_key=TRUE)
  combined_results_with_counts_long=merge(combined_results_with_counts_long,colData(dds_new)[,c("sample","condition")])
  
  
  
  
  pca_plots <- list()
  
  # Create a loop to generate PCA plots for ntop values from 10 to 100 by 10
  for (ntop in seq(10, 100, by = 10)) {
    # Perform variance-stabilizing transformation
    rld <- vst(dds_new, blind = FALSE,fitType ="local")
    
    # Generate PCA data
    pcaData_raw <- plotPCA(rld, intgroup=c("condition", "replicate"), returnData=TRUE, ntop=ntop)
    percentVar_raw <- round(100 * attr(pcaData_raw, "percentVar"))
    
    # Create the PCA plot and add labels
    pca_plot_raw <- ggplot(pcaData_raw, aes(x = PC1, y = PC2, color = condition, shape = replicate)) +
      geom_point(size = 3) +
      geom_text_repel(aes(label = replicate), vjust = 2, size = 3) +
      theme_minimal() +
      labs(title = paste("PCA plot for DESeq analysis. Top ", ntop, " genes. \nSet:",subset_option, sep=""),
           x = paste0("PC1: ", percentVar_raw[1], "% variance"),
           y = paste0("PC2: ", percentVar_raw[2], "% variance"))
    
    # Save the plot to the list
    pca_plots[[paste0("ntop_", ntop)]] <- pca_plot_raw
  }
  #?ggarrange
  pca_plot_corr10_100=ggarrange(plotlist = pca_plots, ncol = 2, nrow = 5,common.legend = TRUE)
  
  
  
  pca_plots <- list()
  
  # Create a loop to generate PCA plots for ntop values from 10 to 100 by 10
  for (ntop in seq(100, 1000, by = 100)) {
    # Perform variance-stabilizing transformation
    rld <- vst(dds_new, blind = FALSE,fitType ="local")
    
    # Generate PCA data
    pcaData_raw <- plotPCA(rld, intgroup=c("condition", "replicate"), returnData=TRUE, ntop=ntop)
    percentVar_raw <- round(100 * attr(pcaData_raw, "percentVar"))
    
    # Create the PCA plot and add labels
    pca_plot_raw <- ggplot(pcaData_raw, aes(x = PC1, y = PC2, color = condition, shape = replicate)) +
      geom_point(size = 3) +
      geom_text_repel(aes(label = replicate), vjust = 2, size = 3) +
      theme_minimal() +
      labs(title = paste("PCA plot for DESeq analysis. Top ", ntop, " genes. \nSet:",subset_option, sep=""),
           x = paste0("PC1: ", percentVar_raw[1], "% variance"),
           y = paste0("PC2: ", percentVar_raw[2], "% variance"))
    
    # Save the plot to the list
    pca_plots[[paste0("ntop_", ntop)]] <- pca_plot_raw
  }
  #?ggarrange
  pca_plot_corr100_1000=ggarrange(plotlist = pca_plots, ncol = 2, nrow = 5,common.legend = TRUE)
  
  
  return(list(dds=dds_new,
              PCAplot_data=vsd,
              pvalue_good=pvalue_good,
              combined_results_with_norm_counts=combined_results_with_norm_counts,
              PCAggplot10_100=pca_plot_corr10_100,
              PCAggplot100_1000=pca_plot_corr100_1000,
              DEGresults = combined_results,
              filteredDEG_all=test_DGE,
              filteredDEG_final=final_test_DGE_no_threshold,
              sigOE=sigOE_long,
              all_genes=combined_results_with_counts_long,
              # up_dpa3_venn_not_filtered=up_dpa3_venn_not_filtered,
              # up_dpa5_venn_not_filtered=up_dpa5_venn_not_filtered,
              # down_dpa3_venn_not_filtered=down_dpa3_venn_not_filtered,
              # down_dpa5_venn_not_filtered=down_dpa5_venn_not_filtered,
              up_dpa3_venn_filtered=up_dpa3_venn_filtered,
              up_dpa5_venn_filtered=up_dpa5_venn_filtered,
              down_dpa3_venn_filtered=down_dpa3_venn_filtered,
              down_dpa5_venn_filtered=down_dpa5_venn_filtered,
              up_dpa3=up_dpa3,
              up_dpa5=up_dpa5,
              down_dpa3=down_dpa3,
              down_dpa5=down_dpa5
              
  ))
}

# LFC <- 0.58
# Fold_Change <- 2^LFC
# Fold_Change #1,49
#Check following:
#dpa3_reduced_Elac21_WT53_GFP23, dpa3_reduced_Elac25_WT53_GFP23
#dpa5_reduced_Elac15_WT52_GFP23, dpa5_reduced_Elac15_WT52_GFP24

#??Deseq
local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23=my_DESeq(counts(all_stringtie_dds_dpa3), subset_option="dpa3_reduced_Elac21_WT53_GFP23",test = "Wald", fitType = "local", reduced = NULL, 
                                                                   minReplicatesForReplace = 7, sfType = "ratio", betaPrior = FALSE,
                                                                   cooksCutoff = TRUE, independentFiltering = TRUE, alpha = 0.1,
                                                                   pAdjustMethod = "BH", parallel = FALSE, minmu = 0.5)


local_Wald_DESeq_corrected_dpa5_reduced_Elac15_WT52_GFP24=my_DESeq(counts(all_stringtie_dds_dpa5), subset_option="dpa5_reduced_Elac15_WT52_GFP24",test = "Wald", fitType = "local", reduced = NULL, 
                                                                   minReplicatesForReplace = 7, sfType = "ratio", betaPrior = FALSE,
                                                                   cooksCutoff = TRUE, independentFiltering = TRUE, alpha = 0.1,
                                                                   pAdjustMethod = "BH", parallel = FALSE, minmu = 0.5)



final_dpa3=local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23
final_dpa5=local_Wald_DESeq_corrected_dpa5_reduced_Elac15_WT52_GFP24



dpa3_reduced_Elac21_WT53_GFP23_venn=ggarrange(local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23$up_dpa3_venn_filtered,
                                              local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23$down_dpa3_venn_filtered)


local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23$down_dpa3$Elac_vs_WT_dpa3[
  local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23$down_dpa3$Elac_vs_WT_dpa3 %in% local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23$down_dpa3$Elac_vs_GFP_dpa3 &
    local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23$down_dpa3$Elac_vs_WT_dpa3 %in% local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23$down_dpa3$GFP_vs_WT_dpa3
]

local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23$filteredDEG_all[local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23$filteredDEG_all$gene %in% 
                                                                            c("SMESG000028264.1","SMESG000035328.1","SMESG000035335.1"),]

local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23$filteredDEG_all

#?annotate_figure
annotate_figure(dpa3_reduced_Elac21_WT53_GFP23_venn,
                top = text_grob("dpa3_reduced_Elac21_WT53_GFP23 (abs(LFC)>0.58)", color = "red", face = "bold", size = 14),
                left = text_grob("Upregulated genes", color = "black",
                                 face = "italic", size = 10),
                right  = text_grob("Downregulated genes", color = "black",
                                   face = "italic", size = 10),
)


dpa5_reduced_Elac15_WT52_GFP24_venn=ggarrange(local_Wald_DESeq_corrected_dpa5_reduced_Elac15_WT52_GFP24$up_dpa5_venn_filtered,
                                              local_Wald_DESeq_corrected_dpa5_reduced_Elac15_WT52_GFP24$down_dpa5_venn_filtered)
#?annotate_figure
annotate_figure(dpa5_reduced_Elac15_WT52_GFP24_venn,
                top = text_grob("dpa5_reduced_Elac15_WT52_GFP24 (abs(LFC)>0.58)", color = "red", face = "bold", size = 14),
                left = text_grob("Upregulated genes", color = "black",
                                 face = "italic", size = 10),
                right  = text_grob("Downregulated genes", color = "black",
                                   face = "italic", size = 10),
)


# ###############################################################################
#________________________________ANNOTATION______________________________________
#################################################################################


final_dpa3=local_Wald_DESeq_corrected_dpa3_reduced_Elac21_WT53_GFP23
final_dpa5=local_Wald_DESeq_corrected_dpa5_reduced_Elac15_WT52_GFP24
# load("E:/Stringtie_anno/SM_anno/final/replicate4_5/anno")
# anno
load_RData_from_github('pfam_swiss_ncbi_merged_only_genes_dedup.RData')
pfam_swiss_ncbi_merged_only_genes_dedup

filteredDEG_final_dpa3=final_dpa3$filteredDEG_final
filteredDEG_final_dpa5=final_dpa5$filteredDEG_final

filteredDEG_final_dpa3=filteredDEG_final_dpa3[filteredDEG_final_dpa3$set=="Elac_vs_WT_dpa3",]
filteredDEG_final_dpa5=filteredDEG_final_dpa5[filteredDEG_final_dpa5$set=="Elac_vs_WT_dpa5",]

final_test_DGE=rbind(filteredDEG_final_dpa3[,c(1:8)],filteredDEG_final_dpa5[,c(1:8)])
nrow(final_test_DGE)#129
#final_test_DGE[final_test_DGE$gene=="SMESG0000037458.1",] #no ELAC
load_RData_from_github('list_GO_gene_correct.RData') #list_GO
list_GO_gene_correct
#remove ELAC from results: SMEST037458003.1, SMEST037458001.1 and "MSTRG.12918.1
#GO:0051082
final_test_DGE

colnames(pfam_swiss_ncbi_merged_only_genes_dedup)
colnames(pfam_swiss_ncbi_merged_only_genes_dedup)=c("gene","Uniprot_protein_name","PFAM_domain_name","NCBI_ID","all_anno","NCBI_only" ,"Annotation")
pfam_swiss_ncbi_merged_only_genes_dedup$all_anno=ifelse(pfam_swiss_ncbi_merged_only_genes_dedup$NCBI_only=="only_NCBI_anno",
                                                        paste(pfam_swiss_ncbi_merged_only_genes_dedup$Annotation,pfam_swiss_ncbi_merged_only_genes_dedup$all_anno, sep="; "),
                                                        pfam_swiss_ncbi_merged_only_genes_dedup$all_anno)
anno=pfam_swiss_ncbi_merged_only_genes_dedup[,c("gene","all_anno")]

head(final_test_DGE)
head(list_GO_gene_correct)
head(names(list_GO_gene_correct))
head(pfam_swiss_ncbi_merged_only_genes_dedup)
length(list_GO_gene_correct) #18522
length(names(list_GO_gene_correct)) #18522
nrow(pfam_swiss_ncbi_merged_only_genes_dedup) #25265
length(na.omit(pfam_swiss_ncbi_merged_only_genes_dedup$Annotation)) #24101


filteredDEG_final_dpa3_annotated <- filteredDEG_final_dpa3[,c(1:17)] %>% 
  dplyr::left_join(anno, by = "gene") %>% 
  dplyr::mutate(gene = paste0(gene, " (", all_anno, ")"))

filteredDEG_final_dpa5_annotated <- filteredDEG_final_dpa5[,c(1:17)] %>% 
  dplyr::left_join(anno, by = "gene") %>% 
  dplyr::mutate(gene = paste0(gene, " (", all_anno, ")"))

# Convert the list to a dataframe
list_GO_gene_correct_df <- do.call(rbind, lapply(names(list_GO_gene_correct), function(gene) {
  data.frame(gene = gene, GO_ID = list_GO_gene_correct[[gene]])
}))

# Reset the row names for clarity
rownames(list_GO_gene_correct_df) <- NULL

# Display the resulting dataframe
print(list_GO_gene_correct_df)

(GO_parent_child=get_parent_nodes(list_GO_gene_correct_df$GO_ID))
GO_parent_child$GO_ID=GO_parent_child$parent_go_id
list_GO_gene_correct_df_with_anno=merge(list_GO_gene_correct_df,GO_parent_child[c("GO_ID","parent_name")],by="GO_ID")
list_GO_gene_correct_df_with_anno=unique(list_GO_gene_correct_df_with_anno)
# Group by gene and combine parent_name values
list_GO_gene_correct_df_with_anno_merged <- list_GO_gene_correct_df_with_anno %>%
  group_by(gene) %>%
  summarize(
    #GO_ID = first(GO_ID), # Optionally keep one GO_ID or use unique for combined GO_IDs
    parent_name = paste(unique(parent_name), collapse = "; ")
  )
load_RData_from_github('only_NCBI_DEG_split.RData')
only_NCBI_DEG_split
only_NCBI_DEG_split_merged <- only_NCBI_DEG_split[,c(1,3)] %>%
  group_by(gene) %>%
  summarize(
    #GO_ID = first(GO_ID), # Optionally keep one GO_ID or use unique for combined GO_IDs
    NCBI_ID_explained = paste(unique(NCBI_ID_full), collapse = "; ")
  )


genes_with_description_and_GO=merge(anno,list_GO_gene_correct_df_with_anno_merged,by="gene",all=TRUE)
genes_with_description_and_GO=merge(genes_with_description_and_GO,only_NCBI_DEG_split_merged,by="gene",all=TRUE)
#save(genes_with_description_and_GO,file="E:/Stringtie_anno/SM_anno/final/final_final/genes_with_description_and_GO.RData")
#write.xlsx(genes_with_description_and_GO,file="E:/Stringtie_anno/SM_anno/final/final_final/genes_with_description_and_GO.xlsx")
anno

filteredDEG_final_dpa3_annotated <- filteredDEG_final_dpa3[,c(1:17)] %>% 
  dplyr::left_join(genes_with_description_and_GO, by = "gene")

filteredDEG_final_dpa5_annotated <- filteredDEG_final_dpa5[,c(1:17)] %>% 
  dplyr::left_join(genes_with_description_and_GO, by = "gene") 


#SMESG000059853.1 is 	heat shock protein beta-1
nrow(anno[is.na(anno$Annotation),]) #1163
anno[anno$gene=="SMESG000059853.1",]


# Load the openxlsx package
library(openxlsx)
# Create a new workbook
wb <- createWorkbook()
# Add the first data frame to the workbook
addWorksheet(wb, "dpa3")
#colnames(final_test_DGE_annotated)
writeData(wb, "dpa3", filteredDEG_final_dpa3_annotated)
# Add the second data frame to the workbook
addWorksheet(wb, "dpa5")
writeData(wb, "dpa5",filteredDEG_final_dpa5_annotated )
DEG_genes=unique(c(filteredDEG_final_dpa3_annotated$gene,filteredDEG_final_dpa5_annotated$gene))

pfam_swiss_ncbi_merged_only_genes_dedup[pfam_swiss_ncbi_merged_only_genes_dedup$gene %in% DEG_genes,]

#new_SMEST_GO_list
DGE_set_with_anno=final_test_DGE
#DGE_set_with_anno=read.csv("G:/silencing_with_batch/new_DGE_set_with_annotation.csv", header = TRUE)
head(DGE_set_with_anno)
str(head(list_GO_gene_correct))
typeof(list_GO_gene_correct)
set=unique(DGE_set_with_anno$set)
# [1] "PARN_dpa3_vs_dpa5" "PARN_vs_GFP_dpa3"  "GFP_vs_WT_dpa5"    "GFP_dpa3_vs_dpa5"  "GFP_vs_WT_dpa3"    "PARN_vs_WT_dpa5"   "Elac_vs_WT_dpa3"   "PARN_vs_GFP_dpa5"  "Elac_vs_GFP_dpa5"  "PARN_vs_WT_dpa3"   "Elac_dpa3_vs_dpa5"
# [12] "Elac_vs_WT_dpa5"   "WT_dpa3_vs_dpa5"   "Elac_vs_GFP_dpa3" 
#DGE_set_with_anno$Transcript[DGE_set_with_anno$set=="WT_vs_PARN_dpa5"&DGE_set_with_anno$log2FoldChange>0]
onts = c( "MF", "BP", "CC" )
GO2geneID <- inverseList(list_GO_gene_correct)
GO2geneID$`GO:0051082`
str(head(GO2geneID))
geneNames <- names(list_GO_gene_correct)
head(geneNames)
#setwd("E:/Illumina/PARN_ELAC2_silencing/mRNA/plots/GO/genes_new_new_new")

# Function to ensure an element is a data frame and not empty
ensureDataFrame <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0) {
    # Create an empty data frame with the same structure as your expected output
    return(data.frame(GO.ID = character(), Expressed_genes = character(), stringsAsFactors = FALSE))
  }
  return(x)
}
#for onthology BP
tab = as.list(onts)
names(tab) = onts
for (i in 1:length(set)){
  tab_over = as.list(onts)
  names(tab_over) = onts
  tab_under = as.list(onts)
  names(tab_under) = onts
  for(j in 1:3){
    print(paste("Gene ontology for",set[i],onts[j], sep=" "))
    genes_over=DGE_set_with_anno$gene[DGE_set_with_anno$set==set[i]&DGE_set_with_anno$log2FoldChange>0]
    
    #over expressed
    #(myInterestingGenes <- genes_over)
    (geneList_over <- factor(as.integer(geneNames %in% genes_over)))
    names(geneList_over) <- geneNames
    #str(geneList)
    # Check if geneList_over has only one level
    if (length(levels(geneList_over)) == 1) {
      next # Skip to the next iteration
    }
    GOdata_over <- new("topGOdata", ontology = onts[j], allGenes = geneList_over,
                       annot = annFUN.gene2GO, gene2GO = list_GO_gene_correct)
    #GOdata@ontology
    resultTopGO.elim_over <- runTest(GOdata_over, algorithm = "elim", statistic = "Fisher" )
    resultTopGO.classic_over <- runTest(GOdata_over, algorithm = "classic", statistic = "Fisher" )
    
    gentable_GO_over=GenTable( GOdata_over, Fisher.elim = resultTopGO.elim_over, 
                               Fisher.classic = resultTopGO.classic_over,
                               orderBy = "Fisher.elim" , topNodes = 200)
    
    gentable_GO_filtered_over=gentable_GO_over[as.numeric(gentable_GO_over$Fisher.elim)<0.05,]
    allGO_over = genesInTerm(GOdata_over) ##get all GOs and their genes from the topgo result.
    #extract just expressed genes
    ANOTATION_over = lapply(allGO_over,function(x) x[x %in% genes_over] )
    #merge elements
    new_ANOTATION_over=lapply(seq(ANOTATION_over[lapply(ANOTATION_over,length)>0]),function(i) paste(ANOTATION_over[lapply(ANOTATION_over,length)>0][[i]], collapse = " "))
    names(new_ANOTATION_over)=names(ANOTATION_over[lapply(ANOTATION_over,length)>0])
    if (nrow(gentable_GO_filtered_over) > 0) {
      df_ANOTATION_over=ldply (new_ANOTATION_over, data.frame)
      colnames(df_ANOTATION_over)=c("GO.ID","Expressed_genes")
      
      # Add genes into gentable
      gentable_GO_merged_over=merge(gentable_GO_filtered_over,df_ANOTATION_over,by="GO.ID",all.x=TRUE)
      gentable_GO_merged_over$sets=set[i]
      gentable_GO_merged_over$GO_class=onts[j]
      gentable_GO_merged_over$expression="overexpressed"
    } else {
      # Handle the case where there are no rows (e.g., create an empty data frame with the same structure)
      gentable_GO_merged_over=data.frame(GO.ID = "NA", Expressed_genes = "NA", 
                                         sets = set[i], GO_class = onts[j], expression = "overexpressed", 
                                         stringsAsFactors = FALSE)
    }
    tab_over[[j]]=gentable_GO_merged_over
    # printGraph(GOdata_over, resultTopGO.elim_over, firstSigNodes =ifelse(nrow(gentable_GO_filtered_over)<10, nrow(gentable_GO_filtered_over), 10),
    #            fn.prefix = paste(set[i],onts[j],"enrichment_overexpressed",sep="_"), useInfo = "all", pdfSW = TRUE)
    # printGraph(GOdata_over, resultTopGO.elim_over, firstSigNodes =ifelse(nrow(gentable_GO_filtered_over)<20, nrow(gentable_GO_filtered_over), 20),
    #            fn.prefix = paste(set[i],onts[j],"enrichment_overexpressed",sep="_"), useInfo = "all", pdfSW = TRUE)
    
    
    #under expressed
    #(myInterestingGenes <- genes_under)
    genes_under=DGE_set_with_anno$gene[DGE_set_with_anno$set==set[i]&DGE_set_with_anno$log2FoldChange<0]
    (geneList_under <- factor(as.integer(geneNames %in% genes_under)))
    names(geneList_under) <- geneNames
    #str(geneList)
    # Check if geneList_over has only one level
    if (length(levels(geneList_under)) == 1) {
      next # Skip to the next iteration
    }
    GOdata_under <- new("topGOdata", ontology = onts[j], allGenes = geneList_under,
                        annot = annFUN.gene2GO, gene2GO = list_GO_gene_correct)
    #GOdata@ontology
    resultTopGO.elim_under <- runTest(GOdata_under, algorithm = "elim", statistic = "Fisher" )
    resultTopGO.classic_under <- runTest(GOdata_under, algorithm = "classic", statistic = "Fisher" )
    gentable_GO_under=GenTable( GOdata_under, Fisher.elim = resultTopGO.elim_under, 
                                Fisher.classic = resultTopGO.classic_under,
                                orderBy = "Fisher.elim" , topNodes = 200)
    
    gentable_GO_filtered_under=gentable_GO_under[as.numeric(gentable_GO_under$Fisher.elim)<0.05,]
    allGO_under = genesInTerm(GOdata_under) ##get all GOs and their genes from the topgo result.
    #extract just expressed genes
    ANOTATION_under = lapply(allGO_under,function(x) x[x %in% genes_under] )
    #merge elements
    new_ANOTATION_under=lapply(seq(ANOTATION_under[lapply(ANOTATION_under,length)>0]),function(i) paste(ANOTATION_under[lapply(ANOTATION_under,length)>0][[i]], collapse = " "))
    names(new_ANOTATION_under)=names(ANOTATION_under[lapply(ANOTATION_under,length)>0])
    if (nrow(gentable_GO_filtered_under) > 0) {
      df_ANOTATION_under=ldply (new_ANOTATION_under, data.frame)
      colnames(df_ANOTATION_under)=c("GO.ID","Expressed_genes")
      
      # Add genes into gentable
      gentable_GO_merged_under=merge(gentable_GO_filtered_under,df_ANOTATION_under,by="GO.ID",all.x=TRUE)
      gentable_GO_merged_under$sets=set[i]
      gentable_GO_merged_under$GO_class=onts[j]
      gentable_GO_merged_under$expression="underexpressed"
    } else {
      # Handle the case where there are no rows
      gentable_GO_merged_under=data.frame(GO.ID = "NA", Expressed_genes = "NA", 
                                          sets = set[i], GO_class = onts[j], expression = "underexpressed", 
                                          stringsAsFactors = FALSE)
    }
    tab_under[[j]]=gentable_GO_merged_under
    # printGraph(GOdata_under, resultTopGO.elim_under, firstSigNodes =ifelse(nrow(gentable_GO_filtered_under)<10, nrow(gentable_GO_filtered_under), 10),
    #            fn.prefix = paste(set[i],onts[j],"enrichment_underexpressed",sep="_"), useInfo = "all", pdfSW = TRUE)
    # printGraph(GOdata_under, resultTopGO.elim_under, firstSigNodes =ifelse(nrow(gentable_GO_filtered_under)<20, nrow(gentable_GO_filtered_under), 20),
    #            fn.prefix = paste(set[i],onts[j],"enrichment_underexpressed",sep="_"), useInfo = "all", pdfSW = TRUE)
    
    
    
  }
  tab_over <- lapply(tab_over, ensureDataFrame)
  tab_under <- lapply(tab_under, ensureDataFrame)
  
  topGOResults_over <- plyr::rbind.fill(tab_over)
  topGOResults_under <- plyr::rbind.fill(tab_under)
  
  # (topGOResults_over <- rbind.fill(tab_over))
  # (topGOResults_under <- rbind.fill(tab_under))
  
  assign(paste(set[i],"topGOResults_over",sep="_"),topGOResults_over)
  assign(paste(set[i],"topGOResults_under",sep="_"),topGOResults_under)
}
################################################################################

all_GO=rbind(
  Elac_vs_WT_dpa3_topGOResults_over,Elac_vs_WT_dpa3_topGOResults_under,
  Elac_vs_WT_dpa5_topGOResults_over,Elac_vs_WT_dpa5_topGOResults_under
  
)
all_GO=na.omit(all_GO)
#sigGenes(GOdata_under)

(GO_parent_child=get_parent_nodes(all_GO$GO.ID))
GO_parent_child$GO.ID=GO_parent_child$parent_go_id
all_GO_elim=merge(all_GO,GO_parent_child[c("GO.ID","parent_name")],by="GO.ID")
all_GO_elim=unique(all_GO_elim)
all_GO_elim=all_GO_elim[-2]
#rename(all_GO_elim, GO_term = parent_name)
names(all_GO_elim)[names(all_GO_elim) == "parent_name"] <- "GO_term"
all_GO_elim$GO_term=sub(",", ".", all_GO_elim$GO_term)
all_GO_elim$GO_term=sub(",", ".", all_GO_elim$GO_term)
all_GO_elim$GO_term=sub(",", ".", all_GO_elim$GO_term)
all_GO_elim$GO_term=sub(",", ".", all_GO_elim$GO_term)
all_GO_elim$GO_term=sub(",", ".", all_GO_elim$GO_term)
#Check 3`-tRNA GO:
all_GO_elim[grep("tRNA",all_GO_elim$GO_term),] #SMESG000037458.1
#all_GO_elim=all_GO_elim[all_GO_elim$Expressed_genes != "SMESG000037458.1",] #no ELAC anyway
all_GO_elim_unchanged=all_GO_elim

all_GO_elim$Fisher.elim_log=log(as.numeric(all_GO_elim$Fisher.elim))
unique(all_GO_elim$sets)
all_GO_elim$dpa=substr(all_GO_elim$sets,nchar(all_GO_elim$sets)-3,nchar(all_GO_elim$sets))
ELAC_GO_elim=all_GO_elim
ELAC_GO_elim_vs_WT=ELAC_GO_elim[ELAC_GO_elim$sets %in% c("Elac_vs_WT_dpa3","Elac_vs_WT_dpa5"),]
#ELAC_GO_elim_vs_GFP=ELAC_GO_elim[ELAC_GO_elim$sets %in% c("Elac_vs_GFP_dpa3","Elac_vs_GFP_dpa5"),]
#Check spindel
ELAC_GO_elim_vs_WT[str_detect(ELAC_GO_elim_vs_WT$GO_term, regex("spindle", ignore_case = TRUE)),]
ELAC_GO_elim_vs_WT[str_detect(ELAC_GO_elim_vs_WT$GO_term, regex("mito", ignore_case = TRUE)),]
ELAC_GO_elim_vs_WT[str_detect(ELAC_GO_elim_vs_WT$GO_term, regex("microtu", ignore_case = TRUE)),]
ELAC_GO_elim_vs_WT[str_detect(ELAC_GO_elim_vs_WT$GO_term, regex("chrom", ignore_case = TRUE)),]
ELAC_GO_elim_vs_WT[str_detect(ELAC_GO_elim_vs_WT$GO_term, regex("cycle", ignore_case = TRUE)),]
ELAC_GO_elim_vs_WT[str_detect(ELAC_GO_elim_vs_WT$GO_term, regex("mios", ignore_case = TRUE)),]
ELAC_GO_elim_vs_WT[str_detect(ELAC_GO_elim_vs_WT$GO_term, regex("divis", ignore_case = TRUE)),]
ELAC_GO_elim_vs_WT[str_detect(ELAC_GO_elim_vs_WT$GO_term, regex("phase", ignore_case = TRUE)),]
ELAC_GO_elim_vs_WT[str_detect(ELAC_GO_elim_vs_WT$GO_term, regex("synap", ignore_case = TRUE)),]
ELAC_GO_elim_vs_WT[str_detect(ELAC_GO_elim_vs_WT$GO_term, regex("neur", ignore_case = TRUE)),]

ELAC_GO_elim_vs_WT$GO_term_len=nchar(ELAC_GO_elim_vs_WT$GO_term)

ELAC_GO_elim_vs_WT[ELAC_GO_elim_vs_WT$GO_class=="MF",]
ELAC_GO_elim_vs_WT[ELAC_GO_elim_vs_WT$GO_term_len>80,]


ELAC_GO_elim_vs_WT$GO_term=chartr(".",",",ELAC_GO_elim_vs_WT$GO_term)
ELAC_GO_elim_vs_WT$GO_term <- str_wrap(ELAC_GO_elim_vs_WT$GO_term, width = 60)
library(openxlsx)
# write.xlsx(ELAC_GO_elim_vs_WT,
#            "E:/Stringtie_anno/SM_anno/final/replicate4_5/ELAC_GO_elim_vs_WT_new.xlsx")


################################################################################
#Annotated genes
ELAC_GO_elim_vs_WT=read.xlsx( "E:/Stringtie_anno/SM_anno/final/replicate4_5/ELAC_GO_elim_vs_WT_new.xlsx")
head(ELAC_GO_elim_vs_WT)
anno
#Convert the anno data frame to a named vector for easier lookup
anno_vector <- setNames(anno$Annotation, anno$gene)

#Function to annotate genes
annotate_genes <- function(genes_str) {
  genes <- strsplit(genes_str, " ")[[1]]
  annotated <- sapply(genes, function(gene) {
    annotation <- anno_vector[gene]
    if (!is.na(annotation)) {
      paste0(gene, " (", annotation, ")")
    } else {
      gene
    }
  })
  paste(annotated, collapse = " ")
}

# Add the annotated_genes column to ELAC_GO_elim_vs_WT
ELAC_GO_elim_vs_WT <- ELAC_GO_elim_vs_WT %>%
  mutate(annotated_genes = sapply(Expressed_genes, annotate_genes))
# write.xlsx(ELAC_GO_elim_vs_WT,
#            "E:/Stringtie_anno/SM_anno/final/replicate4_5/ELAC_GO_elim_vs_WT_annotated.xlsx")

filteredDEG_final_dpa3[,c(1:17)]
filteredDEG_final_dpa5[,c(1:17)]
?write.xlsx

#Check MSTRG.27061, SMESG000036376.1, 


pfam_swiss_ncbi_merged_only_genes_dedup[pfam_swiss_ncbi_merged_only_genes_dedup$gene=="MSTRG.27061",]
pfam_swiss_ncbi_merged_only_genes_dedup[pfam_swiss_ncbi_merged_only_genes_dedup$gene=="SMESG000036376.1",]



#___________________________________BP__________________________________________
ELAC_GO_elim_vs_WT=read.xlsx("D:/Elac2/final_results/tables/ELAC_GO_elim_vs_WT_annotated.xlsx")
DEG_in_scRNA=read.xlsx("D:/Elac2/final_results/tables/filteredDEG_final_annotated_with_GO_scnRNA.xlsx")
DEG_in_scRNA$gene=sapply(strsplit(DEG_in_scRNA$gene, " "), `[`, 1)
unique(DEG_in_scRNA$Generalized.cell.population.group)



head(DEG_in_scRNA)

head(DEG_in_scRNA)
head(ELAC_GO_elim_vs_WT)
colnames(DEG_in_scRNA)
DEG_in_scRNA$Expression=ifelse(DEG_in_scRNA$log2FoldChange>0,"Overexpressed","Underexpressed")
# Assign colors to cell types 

custom_colors <- c(
  "Parenchymal Cells"             = scales::alpha("#1b9e77", 0.6),  
  "All cells"                     = scales::alpha("salmon", 0.6), 
  "Secretory and parapharyngeal cells" = scales::alpha("#7570b3", 0.6), 
  "Goblet cells"                  = scales::alpha("#e7298a", 0.6),  
  "Undetermined"                  = scales::alpha("lightgrey", 0.6),  
  "Neural cells"                  = scales::alpha("#e6ab02", 0.6), 
  "Protonephridia cells"          = scales::alpha("#a6761d", 0.6),  
  "Epidermal cells"               = scales::alpha("#666666", 0.6),  
  "Muscle cells"                  = scales::alpha("#1f78b4", 0.6),  
  "Neoblasts"                     = scales::alpha("#6a3d9a", 0.6),  
  "Phagocytes"                    = scales::alpha("#b15928", 0.6),  
  "Secretory cells"               = scales::alpha("blue", 0.6),  
  "Goblet and intestine cells"     = scales::alpha("#b2df8a", 0.6),  
  "Pigment cells"                 = scales::alpha("#fdbf6f", 0.6)   
)


DEG_in_scRNA_table=as.data.frame(table(DEG_in_scRNA$set,DEG_in_scRNA$Generalized.cell.population.group,DEG_in_scRNA$Expression))
colnames(DEG_in_scRNA_table)=c("dpa","Cell_type","Expression","Freq")
DEG_in_scRNA_table$dpa=ifelse(DEG_in_scRNA_table$dpa=="Elac_vs_WT_dpa3","dpa3","dpa5")
DEG_in_scRNA_table=DEG_in_scRNA_table[DEG_in_scRNA_table$Freq>0,]
# Modify the dataset to work with facets
DEG_in_scRNA_table <- DEG_in_scRNA_table %>%
  mutate(Expression_colored = ifelse(Expression == "Overexpressed",
                                     "<span style='color:darkgreen;'>Overexpressed</span>",
                                     "<span style='color:darkred;'>Underexpressed</span>")) %>%
  mutate(facet_label = paste0( Expression_colored,"<br> genes at ", dpa))
head(DEG_in_scRNA_table)
DEG_in_scRNA_table$Cell_type
DEG_in_scRNA_table[DEG_in_scRNA_table$Expression=="Overexpressed" & DEG_in_scRNA_table$dpa=="dpa3",]

# Ensure that Expressed_genes in ELAC_GO_elim_vs_WT is a character vector
ELAC_GO_elim_vs_WT$Expressed_genes <- as.character(ELAC_GO_elim_vs_WT$Expressed_genes)

# Split Expressed_genes in ELAC_GO_elim_vs_WT into individual DEGs 
ELAC_GO_elim_vs_WT$DEGs_list <- strsplit(ELAC_GO_elim_vs_WT$Expressed_genes, " ")


deg_lookup <- setNames(DEG_in_scRNA$Generalized.cell.population.group, DEG_in_scRNA$gene)

# Function to retrieve the cell population group for each DEG in Expressed_genes
get_cell_population_group <- function(degs) {

  groups <- deg_lookup[degs]

  if (length(groups) > 0 && any(!is.na(groups))) {
    return(paste(unique(groups), collapse = ", "))
  } else {
    return(NA)
  }
}
library(ggtext) 
#install.packages("ggtext")
# Apply the function to each list of DEGs in ELAC_GO_elim_vs_WT
ELAC_GO_elim_vs_WT$Generalized.cell.population.group <- sapply(ELAC_GO_elim_vs_WT$DEGs_list, get_cell_population_group)

# Remove the temporary DEGs_list column
ELAC_GO_elim_vs_WT$DEGs_list <- NULL

# View the updated table
head(ELAC_GO_elim_vs_WT)

unique(ELAC_GO_elim_vs_WT$Generalized.cell.population.group)





#exclude plant and fungus related terms
ELAC_GO_elim_vs_WT$GO_term
not_relevant_terms=c("structural constituent of cell wall","fungal-type cell wall","root meristem growth","auxin polar transport",
                     "regulation of auxin polar transport","chlorophyll metabolic process", "chlorophyll biosynthetic process", 
                     "flavonoid metabolic process", "alkaloid metabolic process", "alkaloid catabolic process", "meristem development", 
                     "meristem growth",  "negative regulation of cytokinin-activated signaling pathway",
                     "flavonol 3-sulfotransferase activity","cell wall")
ELAC_GO_elim_vs_WT=ELAC_GO_elim_vs_WT[!(ELAC_GO_elim_vs_WT$GO_term %in% not_relevant_terms),]
unique(ELAC_GO_elim_vs_WT$GO_term)
library(dplyr)

colnames(ELAC_GO_elim_vs_WT)

ELAC_GO_elim_vs_WT=ELAC_GO_elim_vs_WT %>% mutate(cluster = case_when(
  GO_term %in% c("'de novo' posttranslational protein folding", "'de novo' protein folding", "chaperone-mediated protein complex assembly",
                 "chaperone-mediated protein folding", "protein insertion into mitochondrial outer membrane", 
                 "protein maturation by protein folding", "protein palmitoylation", "positive regulation of protein processing",
                 "N-glycan processing to lysosome", "positive regulation of protein maturation",
                 "telomerase holoenzyme complex assembly") ~ "Protein Folding & Processing",
  
  GO_term %in% c("3'-phosphoadenosine 5'-phosphosulfate metabolic process", "aflatoxin metabolic process", "alkaloid catabolic process", 
                 "alkaloid metabolic process", "ammonia assimilation cycle", "androgen metabolic process", "ceramide metabolic process", 
                 "cholesterol biosynthetic process", "cholesterol efflux", "cholesterol homeostasis", "cholesterol metabolic process", 
                 "cholesterol storage", "cobalamin biosynthetic process", "carboxylic acid biosynthetic process", "ethanol catabolic process", 
                 "flavonoid metabolic process", "glutamate catabolic process", "glutamine biosynthetic process", 
                 "glutamine family amino acid biosynthetic process", "glutamine family amino acid catabolic process", 
                 "dicarboxylic acid catabolic process", "fatty acid catabolic process", "monocarboxylic acid biosynthetic process", 
                 "NAD catabolic process", "pyruvate biosynthetic process", "reactive oxygen species metabolic process", 
                 "retinoic acid metabolic process", "retinol metabolic process", "prostaglandin biosynthetic process", 
                 "prostaglandin metabolic process", "prostanoid biosynthetic process", "prostanoid metabolic process", 
                 "protoporphyrinogen IX biosynthetic process", "protoporphyrinogen IX metabolic process", "vitamin A metabolic process", 
                 "vitamin D catabolic process", "xenobiotic catabolic process", "xenobiotic metabolic process", 
                 "doxorubicin metabolic process","glucose metabolic process", "glucose catabolic process", "glucose catabolic process to pyruvate",
                 "canonical glycolysis", "glycolytic process through glucose-6-phosphate", "diacylglycerol biosynthetic process", 
                 "icosanoid metabolic process", "isoprenoid metabolic process", "monoterpenoid metabolic process", "NADH regeneration", 
                 "chlorophyll biosynthetic process", "chlorophyll metabolic process", 
                 "purine ribonucleoside bisphosphate metabolic process",
                 "pyruvate kinase activity", "carbohydrate phosphorylation", "positive regulation of gluconeogenesis", 
                 "glutamate decarboxylase activity", "oligosaccharide catabolic process", "regulation of isoprenoid metabolic process",
                 "sulfation", "positive regulation of organic acid transport", "RNA-dependent DNA biosynthetic process",
                 "1-acyl-2-lysophosphatidylserine acylhydrolase activity", "1-alpha.25-dihydroxyvitamin D3 23-hydroxylase activity", 
                 "acylglycerol O-acyltransferase activity", "aryl sulfotransferase activity", "beta-N-acetylgalactosaminidase activity", 
                 "cobyrinic acid a.c-diamide synthase activity", "DNA N-glycosylase activity", "DNA-(apurinic or apyrimidinic site) endonuclease activity",
                 "DNA polymerase binding", "flavonol 3-sulfotransferase activity", "glucan endo-1.3-alpha-glucosidase activity", 
                 "glucocorticoid receptor binding", "glutamate-ammonia ligase activity", "hydroxymethylbilane synthase activity", 
                 "lipid binding", "lysophospholipase activity", "non-membrane spanning protein tyrosine kinase activity",
                 "non-membrane spanning protein tyrosine phosphatase activity", 
                 "phosphatidylserine 1-acylhydrolase activity", "phospholipase A1 activity", "phospholipase binding", 
                 "phosphoprotein binding", "phosphoric ester hydrolase activity", "phosphotransferase activity. for other substituted phosphate\ngroups", 
                 "phosphatidyl phospholipase B activity", "protein-cysteine S-acyltransferase activity", "protein-cysteine S-palmitoyltransferase activity",
                 "protein phosphatase 5 binding", "pyrimidine nucleotide binding", "quinine 3-monooxygenase activity", 
                 "RNA-directed 5'-3' RNA polymerase activity", "S-acyltransferase activity", "testosterone 16-alpha-hydroxylase activity", 
                 "testosterone 16-beta-hydroxylase activity", "testosterone 6-beta-hydroxylase activity",  
                 "serine-type endopeptidase inhibitor activity", "peptidyl-tyrosine dephosphorylation", "negative regulation of peptidyl-tyrosine phosphorylation", 
                 "positive regulation of phosphoprotein phosphatase activity", "uroporphyrinogen-III synthase activity",
                 "thiolester hydrolase activity", "acid-ammonia (or amide) ligase activity", "ammonia ligase activity",
                 "peptidase inhibitor activity", "negative regulation of peptidase activity", "glutamic-type peptidase activity",
                 "peptidyl-pyrromethane cofactor linkage", "positive regulation of arachidonic acid secretion",
                 "regulation of arachidonic acid secretion",
                 "peptidyl-tyrosine dephosphorylation involved in inactivation of protein kinase activity",
                 "peptidyl-tyrosine dephosphorylation involved in inactivation\nof protein kinase activity",
                 "1.8-cineole 2-exo-monooxygenase activity", "aldo-keto reductase (NADP) activity", "caffeine oxidase activity",
                 "demethylase activity", "ferroxidase activity", "fatty acid omega-hydroxylase activity", "arachidonic acid 11.12-epoxygenase activity",
                 "arachidonic acid 14.15-epoxygenase activity", "arachidonic acid 5.6-epoxygenase activity", "arachidonic acid omega-hydroxylase activity",
                 "hydroperoxy icosatetraenoate isomerase activity", "linoleic acid epoxygenase activity", "long-chain fatty acid omega-hydroxylase activity",
                 "NADPH-hemoprotein reductase activity", "oxidoreductase activity. acting on CH or CH2 groups", 
                 "oxidoreductase activity. acting on CH or CH2 groups. quinone\nor similar compound as acceptor", 
                 "oxidoreductase activity. acting on NAD(P)H", "oxidoreductase activity. acting on NAD(P)H. heme protein as\nacceptor", 
                 "oxidoreductase activity. acting on paired donors. with\nincorporation or reduction of molecular oxygen. NAD(P)H as\none donor. and incorporation of one atom of oxygen", 
                 "oxidoreductase activity. acting on paired donors. with\nincorporation or reduction of molecular oxygen. reduced\niron-sulfur protein as one donor. and incorporation of one\natom of oxygen",
                 "quinine 3-monooxygenase activity", "oxidative demethylation", "prostaglandin-E synthase activity", 
                 "nitric-oxide synthase regulator activity","inositol catabolic process", 
                 "positive regulation of nitric-oxide synthase activity","inositol oxygenase activity",
                 "endopeptidase inhibitor activity", "endopeptidase regulator activity", "phosphatidylinositol 3-kinase regulatory subunit binding",
                 "mitogen-activated protein kinase binding", "glucocorticoid receptor binding", "Rho GDP-dissociation inhibitor binding", 
                 "protein phosphatase 5 binding", "protein phosphorylated amino acid binding", "GDP-dissociation inhibitor binding", 
                 "phosphatidylinositol 3-kinase binding", "phospholipase binding", "phosphotyrosine residue binding", 
                 "RNA polymerase II CTD heptapeptide repeat phosphatase\nactivity", "STAT family protein binding", 
                 "phosphatidylinositol phospholipase activity", "phosphatidylserine phospholipase activity", 
                 "phospholipase binding", "gamma-catenin binding","phosphatidylinositol binding","phosphatidylinositol 3-kinase binding", 
                 "phosphatidylinositol 3-kinase regulatory subunit binding", "phosphatidylinositol phospholipase C activity", 
                 "peptidyl-tyrosine phosphorylation", "protein localization to nucleolus",
                 "collagen trimer", "H4/H2A histone acetyltransferase complex", "HSP90-CDC37 chaperone complex", 
                 "glucosidase II complex", "NuA4 histone acetyltransferase complex", "serine-type endopeptidase complex", 
                 "pyruvate kinase complex", "transcription preinitiation complex", "UDP-N-acetylglucosamine-lysosomal-enzyme\nN-acetylglucosaminephosphotransferase complex",
                 "potassium exchanging ATPase complex", "UDP-N-acetylglucosamine-lysosomal-enzyme\nN-acetylglucosaminephosphotransferase complex", 
                 "sphingolipid activator protein activity",
                 "UDP-N-acetylglucosamine-lysosomal-enzyme\nN-acetylglucosaminephosphotransferase activity",
                 "acylglycerol O-acyltransferase activity", "arachidonic acid 11.12-epoxygenase activity", "arachidonic acid 14.15-epoxygenase activity", 
                 "arachidonic acid 5.6-epoxygenase activity", "arachidonic acid omega-hydroxylase activity", "fatty acid omega-hydroxylase activity",
                 "linoleic acid epoxygenase activity", "long-chain fatty acid omega-hydroxylase activity", "steroid hydroxylase activity", 
                 "sterol esterase activity", "vitamin D 24-hydroxylase activity", "vitamin D 25-hydroxylase activity", "vitamin D3 25-hydroxylase activity", 
                 "1-alpha.25-dihydroxyvitamin D3 23-hydroxylase activity", "icosanoid biosynthetic process", 
                 "positive regulation of icosanoid secretion", "regulation of icosanoid secretion",
                 "positive regulation of arachidonic acid secretion","triglyceride lipase activity") ~ "Metabolic Processes",
  
  GO_term %in% c("cellular lipid biosynthetic process", "lipid hydroxylation", "lipid import into cell", "lipid storage", 
                 "lipoprotein catabolic process", "phosphatidylcholine catabolic process", "phosphatidylethanolamine catabolic process",
                 "phosphatidylglycerol metabolic process", "phosphatidylserine metabolic process", "long-chain fatty acid biosynthetic process",
                 "unsaturated fatty acid biosynthetic process", "unsaturated fatty acid metabolic process", "linoleic acid metabolic process",
                 "omega-hydroxylase P450 pathway", "positive regulation of fatty acid transport", "positive regulation of lipid transport",
                 "negative regulation of lipid localization", "negative regulation of lipid storage",
                 "calcium-independent phospholipase A2 activity",
                 "regulation of fatty acid transport") ~ "Lipid Processes",
  
  GO_term %in% c("adenylate cyclase-inhibiting G protein-coupled acetylcholine\nreceptor signaling pathway", 
                 "ARF protein signal transduction", "G protein-coupled acetylcholine receptor signaling pathway", 
                 "insulin receptor recycling", "interferon-gamma-mediated signaling pathway", 
                 "phospholipase C-activating G protein-coupled acetylcholine\nreceptor signaling pathway", 
                 "Fc-gamma receptor signaling pathway", "Fc-gamma receptor signaling pathway involved in phagocytosis", 
                 "Fc receptor mediated stimulatory signaling pathway", "Fc receptor signaling pathway", 
                 "peroxisome proliferator activated\nreceptor signaling pathway", "regulation of ARF protein signal transduction", 
                 "protein kinase C-activating G protein-coupled receptor\nsignaling pathway", 
                 "negative regulation of antigen receptor-mediated signaling\npathway", 
                 "negative regulation of receptor signaling pathway via\nJAK-STAT", "phospholipase C activity", 
                 "positive regulation of tau-protein kinase activity", "negative regulation of tau-protein kinase activity",
                 "regulation of tau-protein kinase activity", 
                 "regulation of MAP kinase activity", "positive regulation of MAPK cascade", "negative regulation of MAPK cascade",
                 "negative regulation of epidermal growth factor receptor\nsignaling pathway",
                 "negative regulation of receptor signaling pathway via STAT", "positive regulation of signaling receptor activity",
                 "negative regulation of signaling receptor activity",
                 "tyrosine phosphorylation of STAT protein", "negative regulation of tyrosine phosphorylation of STAT\nprotein",
                 "positive regulation of tyrosine phosphorylation", "regulation of tyrosine phosphorylation of STAT protein",
                 "negative regulation of interferon-gamma-mediated signaling\npathway",
                 "negative regulation of interleukin-6-mediated signaling\npathway", 
                 "negative regulation of interleukin-4-mediated signaling\npathway", 
                 "negative regulation of interleukin-2-mediated signaling\npathway", 
                 "negative regulation of positive thymic T cell selection",
                 "calcium-mediated signaling using intracellular calcium source", 
                 "peroxisome proliferator activated receptor signaling pathway", 
                 "regulation of peroxisome proliferator activated receptor signaling pathway", 
                 "negative regulation of peroxisome proliferator activated receptor signaling pathway", 
                 "regulation of hepatocyte growth factor receptor signaling pathway",
                 "positive regulation of gene silencing by miRNA", 
                 "positive regulation of RNA interference", "negative regulation of MAP kinase activity", 
                 "negative regulation of protein tyrosine kinase activity",
                 "calcium-mediated signaling using intracellular calcium\nsource", 
                 "regulation of peroxisome proliferator activated receptor\nsignaling pathway",
                 "negative regulation of peroxisome proliferator activated\nreceptor signaling pathway", 
                 "positive regulation of amino acid transport","positive regulation of protein kinase B signaling",
                 "negative regulation of protein kinase B signaling")~ "Signaling Pathways",
  
  GO_term %in% c("defense response", "defense response to Gram-positive bacterium", "complement activation (classical pathway, lectin pathway)",
                 "immune response-regulating cell surface receptor signaling\npathway involved in phagocytosis", 
                 "macrophage chemotaxis", "macrophage homeostasis", "leukocyte degranulation", 
                 "leukocyte migration involved in immune response", "positive regulation of macrophage chemotaxis",
                 "negative regulation of macrophage colony-stimulating factor\nsignaling pathway", 
                 "negative regulation of macrophage differentiation", "positive regulation of Fc-gamma receptor signaling pathway\ninvolved in phagocytosis",
                 "positive regulation of phagocytosis", "negative regulation of inflammatory response to antigenic\nstimulus", 
                 "negative regulation of tumor necrosis factor-mediated\nsignaling pathway", "response to virus", "complement activation. lectin pathway", 
                 "complement activation. classical pathway", "respiratory burst involved in inflammatory response", 
                 "respiratory burst involved in defense response", "respiratory burst after phagocytosis", "respiratory burst", 
                 "negative regulation of cytokinin-activated signaling pathway", "positive regulation of cytotoxic T cell differentiation",
                 "positive regulation of opsonization", "positive regulation of tumor necrosis factor superfamily\ncytokine production", 
                 "positive regulation of tumor necrosis factor production", "negative regulation of macrophage colony-stimulating factor\nsignaling pathway",
                 "regulation of phagocytosis", "negative regulation of phagocytosis", "positive regulation of phagocytosis",
                 "inflammatory response", "negative regulation of type I interferon-mediated signaling pathway",
                 "negative regulation of type I interferon-mediated signaling\npathway") ~ "Immune Response",
  
  GO_term %in% c("angiogenesis", "blood vessel development", "blood vessel endothelial cell differentiation", "blood vessel morphogenesis", 
                 "cardiac right atrium morphogenesis", "bone marrow development", "embryonic digestive tract morphogenesis",
                 "embryonic forelimb morphogenesis", "embryonic hindlimb morphogenesis", "embryonic neurocranium morphogenesis",
                 "embryonic skeletal joint development", "otic vesicle development", "nephric duct development", 
                 "neural crest cell development", "neural crest cell differentiation", "liver morphogenesis", 
                 "lung development", "spleen development", "stem cell development", "sinoatrial node cell development",
                 "sinoatrial valve development", "hindlimb morphogenesis", "pronephric duct development", "forelimb morphogenesis", 
                 "hatching", "meristem development", "meristem growth", "retina layer formation", "embryonic skeletal joint morphogenesis",
                 "pronephric distal tubule morphogenesis", "pattern specification involved in pronephros development", 
                 "pattern specification involved in kidney development", "proximal/distal pattern formation", 
                 "proximal/distal pattern formation involved in nephron\ndevelopment", 
                 "proximal/distal pattern formation involved in pronephric\nnephron development", "ureter morphogenesis",
                 "neuroblast fate commitment", "subthalamus development", "thalamus development", "pericardium development",
                 "motor behavior","respiratory tube development", "renal system pattern specification", 
                 "inductive cell migration", "pancreatic A cell differentiation", "pancreatic A cell development", 
                 "regulation of branching morphogenesis of a nerve",
                 "heart process", "egg activation", "heart contraction", "regulation of heart contraction", 
                 "positive regulation of cardiac muscle contraction",
                 "positive regulation of skeletal muscle fiber development","notochord morphogenesis", 
                 "otolith morphogenesis", "negative regulation of fertilization") ~ "Developmental Processes",
  
  GO_term %in% c("cell proliferation in bone marrow", "chondrocyte development", "common myeloid progenitor cell proliferation",
                 "fat cell proliferation", "myoblast fusion involved in skeletal muscle regeneration",
                 "neural crest cell differentiation", "osteoblast development", "positive regulation of cell population proliferation",
                 "negative regulation of cell population proliferation", "regulation of cell migration", "positive regulation of cell migration",
                 "regulation of macrophage migration", "positive regulation of macrophage migration", "myeloid cell homeostasis", "neutrophil homeostasis",
                 "negative regulation of cell growth involved in cardiac\nmuscle cell development", "regulation of hematopoietic stem cell proliferation",
                 "skeletal muscle tissue regeneration", "regulation of locomotion","negative regulation of epithelial cell migration",
                 "regulation of brown fat cell differentiation", "positive regulation of brown fat cell differentiation", 
                 "regulation of macrophage chemotaxis", "replicative senescence", "negative regulation of cell migration", 
                 "negative regulation of cell motility", "regulation of cell motility") ~ "Cell Proliferation & Differentiation",
  
  GO_term %in% c("behavioral response to pain", "cellular response to amine stimulus", "cellular response to interleukin-4", 
                 "cellular response to interleukin-6", "cellular response to lead ion", 
                 "response to antibiotic", "response to cobalt ion", "response to dietary excess", "response to gravity", 
                 "response to inactivity", "response to muscle inactivity", "response to rapamycin", "response to vitamin A",
                 "response to interleukin-4",
                 "response to interleukin-6", "detection of biotic stimulus", "detection of cell density", 
                 "negative regulation of chemotaxis", "positive regulation of chemotaxis",
                 "negative regulation of inflammatory response", "regulation of granulocyte chemotaxis", "positive regulation of granulocyte chemotaxis",
                 "regulation of cell population proliferation", 
                 "positive regulation of PERK-mediated unfolded protein response", "saliva secretion",
                 "surfactant homeostasis", 
                 "negative regulation of collagen metabolic process", "negative regulation of collagen biosynthetic process",
                 "positive chemotaxis", "negative chemotaxis",  
                 "regulation of fertilization","positive regulation of PERK-mediated unfolded protein\nresponse",
                 "regulation of hepatocyte growth factor receptor signaling\npathway", 
                 "negative regulation of homotypic cell-cell adhesion") ~ "Response to Stimuli & Stress",
  
  GO_term %in% c("blood circulation", "regulation of blood circulation", "vasculature development", "vasculogenesis",
                 "regulation of vascular permeability", "negative regulation of vasculature development", 
                 "negative regulation of blood vessel morphogenesis", "positive regulation of sprouting angiogenesis", 
                 "negative regulation of angiogenesis", "regulation of endothelial cell proliferation", 
                 "negative regulation of endothelial cell proliferation", "regulation of sprouting angiogenesis",
                 "negative regulation of endothelial cell migration", "negative regulation of vascular permeability",
                 "vasculogenesis", "vasculature development", "positive regulation of blood coagulation",
                 "negative regulation of blood coagulation", "regulation of platelet activation", "positive regulation of platelet activation",
                 "negative regulation of platelet activation", "regulation of platelet aggregation", 
                 "negative regulation of platelet aggregation", "positive regulation of coagulation",
                 "positive regulation of hemostasis", "fibrinolysis", "negative regulation of fibrinolysis",
                 "regulation of fibrinolysis", "platelet-derived growth factor receptor binding", 
                 "platelet-derived growth factor receptor signaling pathway",
                 "negative regulation of platelet-derived growth factor receptor signaling pathway",
                 "negative regulation of platelet-derived growth factor receptor-beta signaling pathway", 
                 "low-density lipoprotein particle clearance", "triglyceride-rich lipoprotein particle clearance",
                 "negative regulation of platelet-derived growth factor\nreceptor signaling pathway", 
                 "negative regulation of platelet-derived growth factor\nreceptor-beta signaling pathway") ~ "Blood Related Processes",
  
  GO_term %in% c("T cell apoptotic process", "positive regulation of DNA damage response (signal\ntransduction by p53 class mediator)", 
                 "oxidative stress-induced premature senescence", 
                 "positive regulation of DNA damage response. signal\ntransduction by p53 class mediator",
                 "regulation of endoplasmic reticulum stress-induced intrinsic\napoptotic signaling pathway",
                 "positive regulation of apoptosis", "negative regulation of apoptosis",
                 "positive regulation of endoplasmic reticulum stress-induced\nintrinsic apoptotic signaling pathway",
                 "regulation of endoplasmic reticulum stress-induced intrinsic\napoptotic signaling pathway") ~ "Apoptosis & Cell Death",
  
  GO_term %in% c("androgen metabolic process", "estrogen metabolic process", "hormone metabolic process", "thyroid hormone generation",
                 "thyroid hormone metabolic process", "regulation of hormone levels", "regulation of insulin receptor signaling pathway",
                 "negative regulation of insulin receptor signaling pathway",
                 "negative regulation of cellular response to insulin stimulus", "estrogen 16-alpha-hydroxylase activity",
                 "estrogen 2-hydroxylase activity") ~ "Hormonal & Endocrine Processes",
  
  GO_term %in% c("axon ensheathment", "ensheathment of neurons", "defasciculation of motor neuron axon", "fasciculation of sensory neuron axon",
                 "gamma-aminobutyric acid secretion", "ganglioside catabolic process", "synaptic transmission (GABAergic)", 
                 "synaptic transmission (glutamatergic)", "positive regulation of synaptic transmission (GABAergic, glutamatergic)",
                 "negative regulation of synaptic transmission (GABAergic, glutamatergic)", 
                 "regulation of short-term neuronal synaptic plasticity", "synaptic target inhibition", 
                 "synaptic target recognition", "motor behavior",  "response to interleukin-6", 
                 "negative regulation of synaptic transmission. glutamatergic","neuropeptide signaling pathway", 
                 "positive regulation of synaptic transmission. glutamatergic",
                 "anchoring junction", "cholinergic synapse", "immunological synapse", "protein complex involved in cell-cell adhesion", 
                 "axon terminus", "neuron projection terminus","retinal ganglion cell axon guidance",
                 "kainate selective glutamate receptor complex", "negative regulation of synaptic transmission. GABAergic",
                 "positive regulation of synaptic transmission. GABAergic", "negative regulation of synaptic transmission. glutamatergic", 
                 "positive regulation of synaptic transmission. glutamatergic", "regulation of T cell receptor signaling pathway",
                 "positive regulation of T cell receptor signaling pathway", "negative regulation of T cell receptor signaling pathway", 
                 "positive regulation of antigen receptor-mediated signaling\npathway", "inhibitory postsynaptic potential",
                 "cell adhesion molecule binding", "negative regulation of excitatory synapse assembly", 
                 "regulation of gamma-aminobutyric acid secretion", 
                 "positive regulation of gamma-aminobutyric acid secretion",
                 "myelination") ~ "Neuronal Processes",
  
  GO_term %in% c("acetylcholine receptor activity", "G protein-coupled acetylcholine receptor activity", 
                 "G protein-coupled neurotransmitter receptor activity", "G protein-coupled serotonin receptor activity",
                 "extracellularly glutamate-gated ion channel activity", "kainate selective glutamate receptor activity", 
                 "potassium channel regulator activity", "potassium ion binding", "calcium channel regulator activity", 
                 "ion channel binding", "arsenite transmembrane transporter activity", "acidic amino acid transmembrane transporter activity",
                 "L-glutamate transmembrane transporter activity", "cholesterol transfer activity", "negative regulation of anion channel activity",
                 "positive regulation of anion channel activity", "negative regulation of anion transmembrane transport",
                 "positive regulation of anion transmembrane transport", "positive regulation of protein transport", 
                 "negative regulation of protein transport",
                 "potassium:proton exchanging ATPase complex","positive regulation of anion transport",
                 "glycolipid transport", "intracellular cholesterol transport",
                 "lipid import into cell", "vesicle-mediated transport","intra-Golgi vesicle-mediated transport") ~"Membrane Transport Activities",
  
  GO_term %in% c( "alpha-catenin binding", "beta-catenin binding", "cadherin binding", "calcium-dependent carbohydrate binding", 
                  "cholesterol binding", "delta-catenin binding", "dynein light chain binding",   
                  "fucose binding", "GDP-dissociation inhibitor binding", "glutamate binding", "hormone activity",  
                  "mannose binding", "mitogen-activated protein kinase binding", "nitric-oxide synthase binding", "pyrimidine nucleotide binding",
                  "pyrimidine ribonucleotide binding", "protein phosphorylated amino acid binding", "phosphotyrosine residue binding", 
                  "retinoic acid 4-hydroxylase activity", "STAT family protein binding", "sulfonylurea receptor binding", "Tat protein binding", 
                  "thyroid hormone binding", "vitamin D binding", "transmembrane receptor protein tyrosine phosphatase activity",
                  "UTP binding", "CTP binding", "phospholipid binding",  "phosphoprotein binding", 
                  "TPR domain binding", "oligosaccharide binding", "3'-phosphoadenosine 5'-phosphosulfate binding") ~"Binding Activities",
  
  GO_term %in% c("alkali metal ion binding", "ferric iron binding", "ferrous iron binding", "iron ion binding", "oxygen binding", 
                 "sequestering of iron ion", "sequestering of metal ion","intracellular sequestering of iron ion", 
                 "oxidoreductase activity. acting on metal ions", "oxidoreductase activity. acting on metal ions. oxygen as\nacceptor",
                 "iron binding activity", "oxygen binding activity") ~"Iron, Metal & Oxygen Binding",
  
  GO_term %in% c("1-alpha.25-dihydroxyvitamin D3 23-hydroxylase activity", "anandamide 11.12 epoxidase activity", "anandamide 14.15 epoxidase activity", 
                 "anandamide 8.9 epoxidase activity", "anandamide epoxidase activity", "flavonol 3-sulfotransferase activity", 
                 "thromboxane-A synthase activity", 
                 "urophorphyrinogen-III synthase activity", "toxin activity", "retina layer formation",  
                 "positive regulation of focal adhesion assembly",
                 "negative regulation of focal adhesion assembly", "positive regulation of actin filament polymerization", 
                 "negative regulation of actin filament polymerization", "cell-cell recognition", "prevention of polyspermy", 
                 "defecation", "contact inhibition",
                 "maintenance of CRISPR repeat elements", "envenomation resulting in modulation of process in other organism", 
                 "envenomation resulting in modulation of voltage-gated potassium channel activity in other organism", 
                 "envenomation resulting in modulation of ion channel activity in other organism",
                 "envenomation resulting in negative regulation of voltage-gated potassium channel activity in other organism",
                 "envenomation resulting in modulation of process in other\norganism","snRNA transcription",
                 "envenomation resulting in modulation of voltage-gated\npotassium channel activity in other organism", 
                 "envenomation resulting in modulation of ion channel activity\nin other organism", 
                 "envenomation resulting in negative regulation of\nvoltage-gated potassium channel activity in other organism",
                 "regulation of protein localization to nucleolus", "positive regulation of telomerase activity", 
                 "positive regulation of RNA interference","binding of sperm to zona pellucida", 
                 "negative regulation of binding of sperm to zona pellucida","prevention of polyspermy", 
                 "fusion of sperm to egg plasma membrane involved in single\nfertilization", "sperm-egg recognition"
  ) ~"Miscellaneous Activities",
  
  GO_term %in% c("acrosomal vesicle", "azurophil granule", "azurophil granule lumen", "azurophil granule membrane", "cortical granule", 
                 "extracellular vesicle", "multivesicular body", "platelet alpha granule membrane", "specific granule membrane",
                 "tertiary granule", "tertiary granule membrane",  "endosome lumen", 
                 "platelet alpha granule membrane", "plasma membrane fusion",
                 "multivesicular body", "zymogen granule") ~"Vesicles & Granules",
  
  GO_term %in% c("cell surface", "extracellular membrane-bounded organelle", "membrane microdomain", "membrane raft", 
                 "outer membrane-bounded periplasmic space", "postsynaptic density membrane", 
                 "ruffle membrane",   
                 "cell wall", "primary cell wall", "secondary cell wall") ~"Membrane Structures",
  
  GO_term %in% c("telomere maintenance via telomerase", "telomere maintenance via telomere lengthening",
                 "endocytosis") ~ "Organellar Activities",
  
  GO_term %in% c("extracellular exosome", "extracellular organelle",  "lysosome", 
                 "lysosomal lumen", "primary lysosome", "extracellular region", "extracellular space", 
                 "outer mitochondrial membrane organization", "periplasmic space") ~"Extracellular & Organellar Regions",
  
  GO_term %in% c("dendritic growth cone", "Golgi cis cisterna", "lamellipodium", "rough endoplasmic reticulum", 
                 "sperm connecting piece", "sperm mitochondrial sheath", "nematocyst",
                 "sperm plasma membrane", 
                 "sperm mitochondrial sheath", "sperm connecting piece", "notochord morphogenesis", "otolith formation", 
                 "otolith morphogenesis",
                 "podosome assembly","regulation of podosome assembly","striated muscle myosin thick filament assembly") ~"Specialized Structures",
  TRUE ~ "other"
))




unique(ELAC_GO_elim_vs_WT$cluster)
table(ELAC_GO_elim_vs_WT$cluster)
unique(ELAC_GO_elim_vs_WT$GO_term[ELAC_GO_elim_vs_WT$cluster=="other"])
unique(ELAC_GO_elim_vs_WT$GO_term[ELAC_GO_elim_vs_WT$cluster=="Apoptosis & Cell Death"]) #"T cell apoptotic process"
unique(ELAC_GO_elim_vs_WT$GO_term[ELAC_GO_elim_vs_WT$cluster=="Iron, Metal & Oxygen Binding"]) #"oxygen binding"
write.xlsx(ELAC_GO_elim_vs_WT,file="D:/Elac2/final_results/tables/ELAC_GO_elim_vs_WT_annotated_grouped_and_with_cell_cluster.xlsx")

unique(DEG_in_scRNA$Generalized.cell.population.group)



#WT
str(ELAC_GO_elim_vs_WT)
ELAC_GO_elim_vs_WT_BP=sort(unique(ELAC_GO_elim_vs_WT$GO_term))
ELAC_GO_elim_vs_WT$Significant=as.integer(ELAC_GO_elim_vs_WT$Significant)
length(ELAC_GO_elim_vs_WT_BP) #576
round(length(ELAC_GO_elim_vs_WT_BP)/100) #6 plots
length(ELAC_GO_elim_vs_WT_BP)/6 #96 terms on each plot
head(ELAC_GO_elim_vs_WT)

str(ELAC_GO_elim_vs_WT$Fisher.elim)
str(ELAC_GO_elim_vs_WT$Significant)
ELAC_GO_elim_vs_WT$Fisher.elim=as.numeric(ELAC_GO_elim_vs_WT$Fisher.elim)
#95 terms in each plot

#nTerm=95
#ELAC_GO_elim_vs_WT_BP=ELAC_GO_elim_vs_WT[ELAC_GO_elim_vs_WT$GO_class=="BP",]
table(ELAC_GO_elim_vs_WT$cluster)



unique(DEG_in_scRNA$Generalized.cell.population.group)


# Define your list of cell types
cell_types <- c(
  "Parenchymal Cells", "All cells", "Secretory and parapharyngeal cells", 
  "Goblet cells", "Undetermined", "Neural cells", "Protonephridia cells", 
  "Epidermal cells", "Muscle cells", "Neoblasts", "Phagocytes", 
  "Secretory cells", "Goblet and intestine cells", "Pigment cells"
)

# Loop over each cell type and create a new column
for (cell_type in cell_types) {
  # Create a valid column name by replacing spaces and special characters
  column_name <- make.names(gsub(" ", "_", cell_type))
  
  # Use ifelse to assign values based on matching cell types
  ELAC_GO_elim_vs_WT[[column_name]] <- ifelse(
    grepl(cell_type, ELAC_GO_elim_vs_WT$Generalized.cell.population.group),
    cell_type,
    ""
  )
}


# 
# ELAC_GO_elim_vs_WT$Parenchymal=ifelse(grepl("Parenchymal Cells",
#                                             ELAC_GO_elim_vs_WT$Generalized.cell.population.group),
#                                       "Parenchymal Cells","")
unique(ELAC_GO_elim_vs_WT$cluster)

ELAC_GO_elim_vs_WT$Generalized.cell.population.group=ifelse(grepl("All cells",
                                                                  ELAC_GO_elim_vs_WT$Generalized.cell.population.group),"All cells",
                                                            ELAC_GO_elim_vs_WT$Generalized.cell.population.group)

ELAC_GO_elim_vs_WT$Generalized.cell.population.group=ifelse(ELAC_GO_elim_vs_WT$Generalized.cell.population.group==
                                                              "Goblet and intestine cells, Goblet cells",
                                                            "Goblet and intestine cells",ELAC_GO_elim_vs_WT$Generalized.cell.population.group)
ELAC_GO_elim_vs_WT$Generalized.cell.population.group=ifelse(ELAC_GO_elim_vs_WT$Generalized.cell.population.group %in%
                                                              c("Neural cells, Epidermal cells, Goblet and intestine cells, Goblet cells, Secretory and parapharyngeal cells",
                                                                "Neural cells, Epidermal cells, Secretory and parapharyngeal cells" ),
                                                            "All cells",ELAC_GO_elim_vs_WT$Generalized.cell.population.group)


unique(ELAC_GO_elim_vs_WT$Generalized.cell.population.group)
unique(ELAC_GO_elim_vs_WT$Expressed_genes[ELAC_GO_elim_vs_WT$Generalized.cell.population.group
                                          =="Secretory and parapharyngeal cells, Goblet cells"])


ELAC_GO_elim_vs_WT=ELAC_GO_elim_vs_WT[ELAC_GO_elim_vs_WT$Significant>1,]
# Remove the word "cells" from all values in temp_df
ELAC_GO_elim_vs_WT$Generalized.cell.population.group <- gsub("cells", "", ELAC_GO_elim_vs_WT$Generalized.cell.population.group, ignore.case = TRUE)
ELAC_GO_elim_vs_WT$Generalized.cell.population.group <- gsub("\\s+", " ", ELAC_GO_elim_vs_WT$Generalized.cell.population.group)
ELAC_GO_elim_vs_WT$Generalized.cell.population.group <- trimws(ELAC_GO_elim_vs_WT$Generalized.cell.population.group)
ELAC_GO_elim_vs_WT$Generalized.cell.population.group=ifelse(ELAC_GO_elim_vs_WT$Generalized.cell.population.group=="All",
                                                            "All cells",ELAC_GO_elim_vs_WT$Generalized.cell.population.group)

# Remove any spaces before and after commas
ELAC_GO_elim_vs_WT$Generalized.cell.population.group <- gsub("\\s*,\\s*", ", ", ELAC_GO_elim_vs_WT$Generalized.cell.population.group)



unique(ELAC_GO_elim_vs_WT$Generalized.cell.population.group)
ELAC_GO_elim_vs_WT$Generalized.cell.population.group=ifelse(ELAC_GO_elim_vs_WT$Generalized.cell.population.group=="Goblet and intestine",
                                                            "Goblet, Intestine",ELAC_GO_elim_vs_WT$Generalized.cell.population.group)
ELAC_GO_elim_vs_WT$Generalized.cell.population.group=ifelse(ELAC_GO_elim_vs_WT$Generalized.cell.population.group=="Secretory and parapharyngeal, Goblet, Epidermal",
                                                            "Goblet, Secretory, Parapharyngeal, Epidermal",ELAC_GO_elim_vs_WT$Generalized.cell.population.group)
ELAC_GO_elim_vs_WT$Generalized.cell.population.group=ifelse(ELAC_GO_elim_vs_WT$Generalized.cell.population.group %in% c("Goblet, Secretory and parapharyngeal",
                                                                                                                        "Secretory and parapharyngeal, Goblet"),
                                                            "Goblet, Secretory, Parapharyngeal",ELAC_GO_elim_vs_WT$Generalized.cell.population.group)

ELAC_GO_elim_vs_WT$Generalized.cell.population.group=ifelse(ELAC_GO_elim_vs_WT$Generalized.cell.population.group=="Phagocytes, Neural",
                                                            "Neural, Phagocytes",ELAC_GO_elim_vs_WT$Generalized.cell.population.group)
ELAC_GO_elim_vs_WT$Generalized.cell.population.group=ifelse(ELAC_GO_elim_vs_WT$Generalized.cell.population.group=="Phagocytes, Parenchymal, Neural",
                                                            "Neural, Phagocytes, Parenchymal",ELAC_GO_elim_vs_WT$Generalized.cell.population.group)

ELAC_GO_elim_vs_WT$Generalized.cell.population.group=ifelse(ELAC_GO_elim_vs_WT$Generalized.cell.population.group=="Phagocytes, Parenchymal",
                                                            "Parenchymal, Phagocytes",ELAC_GO_elim_vs_WT$Generalized.cell.population.group)
unique(ELAC_GO_elim_vs_WT$Generalized.cell.population.group)
# Wrap long text 
ELAC_GO_elim_vs_WT$Generalized.cell.population.group <- str_wrap(ELAC_GO_elim_vs_WT$Generalized.cell.population.group, width = 15)

unique(ELAC_GO_elim_vs_WT$cluster)
ELAC_GO_elim_vs_WT$Generalized.cell.population.group
ELAC_GO_elim_vs_WT$cluster=ifelse(ELAC_GO_elim_vs_WT$GO_term %in% c("regulation of hormone levels", 
                                                                    "hormone metabolic process"),"Metabolic Processes",ELAC_GO_elim_vs_WT$cluster)
ELAC_GO_elim_vs_WT$cluster=ifelse(ELAC_GO_elim_vs_WT$cluster %in% c("Membrane Transport Activities", 
                                                                    "Membrane Structures"),"Membrane Structures & Transport",ELAC_GO_elim_vs_WT$cluster)

ELAC_GO_elim_vs_WT$cluster=ifelse(ELAC_GO_elim_vs_WT$cluster=="Vesicles & Granules","Extracellular Regions",ELAC_GO_elim_vs_WT$cluster)
ELAC_GO_elim_vs_WT$cluster=ifelse(ELAC_GO_elim_vs_WT$GO_term %in% c("extracellular region", "extracellular organelle", "extracellular exosome", 
                                                                    "extracellular space"),"Extracellular Regions",ELAC_GO_elim_vs_WT$cluster)
ELAC_GO_elim_vs_WT$cluster=ifelse(ELAC_GO_elim_vs_WT$GO_term %in% c("primary lysosome","lysosomal lumen", "lysosome"),"Organellar Regions & Activities",
                                  ELAC_GO_elim_vs_WT$cluster)
ELAC_GO_elim_vs_WT$cluster=ifelse(ELAC_GO_elim_vs_WT$cluster=="Organellar Activities","Organellar Regions & Activities",ELAC_GO_elim_vs_WT$cluster)


unique(ELAC_GO_elim_vs_WT$cluster)
unique(ELAC_GO_elim_vs_WT$dpa)
unique(ELAC_GO_elim_vs_WT$Generalized.cell.population.group)
unique(ELAC_GO_elim_vs_WT$expression)
unique(ELAC_GO_elim_vs_WT$GO_term)
###############################################################################
#----------------------- Create the pie chart-----------------------------------
custom_colors <- c(
  "Parenchymal cells"             = scales::alpha("#1b9e77", 0.6),  
  "All cells"                     = scales::alpha("salmon", 0.6),  
  "Secretory and parapharyngeal cells" = scales::alpha("#7570b3", 0.6),  
  "Goblet cells"                  = scales::alpha("#e7298a", 0.6),  
  "Undetermined"                  = scales::alpha("lightgrey", 0.6),  
  "Neural cells"                  = scales::alpha("#e6ab02", 0.6),  
  "Protonephridia cells"          = scales::alpha("#a6761d", 0.6),  
  "Epidermal cells"               = scales::alpha("#666666", 0.6),  
  "Muscle cells"                  = scales::alpha("#1f78b4", 0.6),  
  "Neoblasts"                     = scales::alpha("#6a3d9a", 0.6), 
  "Phagocytes"                    = scales::alpha("#b15928", 0.6),  
  "Secretory cells"               = scales::alpha("blue", 0.6),  
  "Goblet and intestine cells"     = scales::alpha("#b2df8a", 0.6),  
  "Pigment cells"                 = scales::alpha("#fdbf6f", 0.6)   
)
DEG_in_scRNA_table$Cell_type <- as.character(DEG_in_scRNA_table$Cell_type)
DEG_in_scRNA_table$Cell_type <- ifelse(DEG_in_scRNA_table$Cell_type=="Parenchymal Cells","Parenchymal cells",DEG_in_scRNA_table$Cell_type)
df_labels <- DEG_in_scRNA_table %>%
  dplyr::group_by(facet_label) %>%             
  dplyr::arrange(Cell_type) %>%               
  dplyr::mutate(
    csum = rev(cumsum(rev(Freq))),      # Cumulative sum 
    pos = Freq / 2 + lead(csum, 1),     # middle of the slice
    pos = if_else(is.na(pos), Freq / 2, pos)  
  ) %>%
  dplyr::ungroup()

#<br>
#?geom_label_repel
head(DEG_in_scRNA_table)
unique(DEG_in_scRNA_table$Cell_type)
(A_fig6=ggplot(DEG_in_scRNA_table, aes(x = "", y = Freq, fill = Cell_type)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = custom_colors) +
    facet_wrap(~ facet_label, nrow = 1, labeller = label_wrap_gen(),scales ="free") +
    theme_void() +
    theme(
      strip.text = element_markdown(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    ) +
    geom_label_repel(
      data = df_labels,
      aes(y = pos, label = Freq),      
      size = 2,
      nudge_x = 1,                                   
      show.legend = FALSE,
      fill = "white",
      color = "black",
      fontface = "bold",
      box.padding = 0.35,
      point.padding = 0.5,
      segment.color = "grey50",
      max.overlaps = Inf                              # Allows all labels to be plotted
    ) +
    labs(
      title = "",
      fill = "Cell type"
    ))+guides(fill=guide_legend(nrow=4,bycol=TRUE))


DEG_in_scRNA_table2 <- DEG_in_scRNA_table %>%
  mutate(
    # wrap at, say, 12 characters per line:
    Cell_type_wrapped = str_wrap(Cell_type, width = 12)
  )

unique(DEG_in_scRNA_table2$Cell_type_wrapped)

custom_colors <- c(
  "Parenchymal\ncells"             = scales::alpha("#1b9e77", 0.6),  
  "All cells"                     = scales::alpha("salmon", 0.6),  
  "Secretory\nand\nparapharyngeal\ncells" = scales::alpha("#7570b3", 0.6), 
  "Goblet cells"                  = scales::alpha("#e7298a", 0.6),  
  "Undetermined"                  = scales::alpha("lightgrey", 0.6),  
  "Neural cells"                  = scales::alpha("#e6ab02", 0.6),  
  "Protonephridia\ncells"          = scales::alpha("#a6761d", 0.6),  
  "Epidermal\ncells"               = scales::alpha("#666666", 0.6),  
  "Muscle cells"                  = scales::alpha("#1f78b4", 0.6),  
  "Neoblasts"                     = scales::alpha("#6a3d9a", 0.6),  
  "Phagocytes"                    = scales::alpha("#b15928", 0.6),  
  "Secretory\ncells"               = scales::alpha("blue", 0.6),  
  "Goblet and\nintestine\ncell"     = scales::alpha("#b2df8a", 0.6),  
  "Pigment\ncells"                 = scales::alpha("#fdbf6f", 0.6)   
)


(A_fig6 <- ggplot(DEG_in_scRNA_table2, 
                  aes(x = "", y = Freq, fill = Cell_type_wrapped)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar(theta = "y") +
    scale_fill_manual(
      values = custom_colors,
      name   = "Cell type"          
    ) +
    facet_wrap(
      ~ facet_label,
      nrow     = 1,
      labeller = label_wrap_gen(),  
      scales   = "free"
    ) +
    theme_void() +
    theme(
      strip.text      = element_markdown(face = "bold", size = 12),
      plot.title      = element_text( face = "bold", size = 12),
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 3, bycol = TRUE)) +
    geom_label_repel(
      data          = df_labels,
      aes(y = pos, label = Freq),
      size          = 2,
      nudge_x       = 1,
      show.legend   = FALSE,
      fill          = "white",
      color         = "black",
      fontface      = "bold",
      box.padding   = 0.35,
      point.padding = 0.5,
      segment.color = "grey50",
      max.overlaps  = Inf
    ))


#?facet_wrap
#_______________________________________________________________________________

#----------------------------Barplots ofr GOs-----------------------------------
head(ELAC_GO_elim_vs_WT)
#ELAC_GO_elim_vs_WT_table=as.data.frame(table(ELAC_GO_elim_vs_WT$cluster))

# Add proportion column
head(ELAC_GO_elim_vs_WT)

ELAC_GO_elim_vs_WT_table <- ELAC_GO_elim_vs_WT %>%
  tidyr::separate_rows(Expressed_genes, sep = " ") %>% # split the space‐separated gene lists into one row per gene
  dplyr::group_by(cluster, dpa, expression) %>%            # Group by cluster, dpa, and expression
  dplyr::summarise(
    unique_GO_IDs = n_distinct(GO.ID),               # Count of unique GO.IDs
    unique_genes    = list(unique(Expressed_genes)), # List of distinct genes
    #sum_Significant = n_distinct(unique_genes), # Sum of Significant genes
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    sum_Significant = lengths(unique_genes)
  )
ELAC_GO_elim_vs_WT_table$expression=ifelse(ELAC_GO_elim_vs_WT_table$expression=="overexpressed","Overexpressed","Underexpressed")
ELAC_GO_elim_vs_WT_table <- ELAC_GO_elim_vs_WT_table %>%
  mutate(Expression_colored = ifelse(expression == "Overexpressed",
                                     "<span style='color:darkgreen;'>Overexpressed</span>",
                                     "<span style='color:darkred;'>Underexpressed</span>")) %>%
  mutate(facet_label = paste0( "<br>",Expression_colored,"<br> genes at ", dpa))

head(ELAC_GO_elim_vs_WT_table)


# Add proportion column
# ELAC_GO_elim_vs_WT_summary_table <- ELAC_GO_elim_vs_WT_table %>%
#   dplyr::group_by(cluster, dpa, expression) %>%
#   dplyr::mutate(
#     proportion = sum_Significant / sum(sum_Significant) * 100
#   ) %>%
#   dplyr::ungroup()
# View the updated table
#print(summary_table)

unique_clusters <- unique(ELAC_GO_elim_vs_WT_table$cluster)
num_clusters <- length(unique_clusters)
custom_colors <- colorRampPalette(brewer.pal(12, "Paired"))(num_clusters)
names(custom_colors) <- unique_clusters
# Calculate label positions for unique_GO_IDs
unique_GO_labels <- ELAC_GO_elim_vs_WT_table %>%
  dplyr::group_by(facet_label) %>%
  dplyr::arrange(cluster) %>%
  dplyr::mutate(
    fraction = unique_GO_IDs / sum(unique_GO_IDs),
    cumulative =rev(cumsum(rev(fraction))),
    pos = cumulative - fraction / 2
  ) %>%
  dplyr::ungroup()

# Calculate label positions for sum_Significant
sum_Significant_labels <- ELAC_GO_elim_vs_WT_table %>%
  dplyr::group_by(facet_label) %>%
  dplyr::arrange(cluster) %>%
  dplyr::mutate(
    fraction = sum_Significant / sum(sum_Significant),
    cumulative = rev(cumsum(rev(fraction))),
    pos = cumulative - fraction / 2
  ) %>%
  dplyr::ungroup()

# Plot Unique GO.IDs
unique_GO_labels[unique_GO_labels$expression=="Overexpressed" & unique_GO_labels$dpa=="dpa5",]
unique_GO_labels$facet_label

# unique_GO_labels <- unique_GO_labels %>%
#   mutate(facet_label = gsub("\n", "<br>", facet_label))


plot_unique_GO <- ggplot(unique_GO_labels, aes(x = 1, y = unique_GO_IDs, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual(values = custom_colors) +
  facet_wrap(~ facet_label, ncol = 2, scales = "free",labeller = label_wrap_gen()) +
  geom_label_repel(
    data = unique_GO_labels,
    aes(x = 1, y = pos, label = unique_GO_IDs),
    inherit.aes = FALSE,
    size = 2,
    nudge_x = 0.3,
    direction = "y",
    segment.color = "grey50",
    show.legend = FALSE,
    box.padding = 0.3,
    point.padding = 0.5,
    max.overlaps = Inf
  ) +
  labs(
    x = "",
    y = "",
    fill = "Clusters",
    title = "Number of GO terms in cluster"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_markdown(face = "bold", size = 12),
    strip.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),       # Remove x-axis text since x is constant
    axis.text.y = element_blank(), 
    #axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 12,hjust = 0.5),
    legend.position = "right",
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
  )


# Plot Sum of Significant Genes
plot_sum_Significant <- ggplot(sum_Significant_labels, aes(x = 1, y = sum_Significant, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual(values = custom_colors) +
  facet_wrap(~ facet_label, ncol = 2, scales = "free",labeller = label_wrap_gen()) +
  geom_label_repel(
    data = sum_Significant_labels,
    aes(x = 1, y = pos, label = sum_Significant),
    inherit.aes = FALSE,
    size = 2,
    nudge_x = 0.3,
    direction = "y",
    segment.color = "grey50",
    show.legend = FALSE,
    box.padding = 0.3,
    point.padding = 0.5,
    max.overlaps = Inf
  ) +
  labs(
    x = "",
    y = "",
    fill = "Clusters",
    title = ""
  ) +
  theme_minimal() +
  theme(
    strip.text = element_markdown(face = "bold", size = 12),
    strip.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 12,hjust = 0.5),
    legend.position = "right",
    plot.margin = unit(c(0,0,0,0), "pt")
  )
library(patchwork)


plot_unique_GO 
# Combine both plots side by side
(B_fig6 <- plot_sum_Significant + theme(plot.margin = margin(0,0,0,0, "cm")))

# ggarrange(A_fig6+guides(fill = guide_legend(nrow=5,bycol=TRUE)),
#           B_fig6,nrow=2, heights = c(0.5,1))+theme(plot.margin = margin(0,0,0,0, "cm")) 

ggarrange(A_fig6,
          B_fig6,nrow=2, heights = c(0.5,1))+theme(plot.margin = margin(0,0,0,0, "cm")) 
#labels=c("A","B"),

print(sum_Significant_labels, n=50)
print(unique_GO_labels, n=50)





head(DEG_in_scRNA)
head(ELAC_GO_elim_vs_WT)
head(ELAC_GO_elim_vs_WT_table)


# gene → GO terms lookup
gene2GO <- ELAC_GO_elim_vs_WT %>%
  tidyr::separate_rows(Expressed_genes, sep = " ") %>%     
  dplyr::rename(                                           
    gene = Expressed_genes
  ) %>%
  dplyr::group_by(                                         
    gene
  ) %>%
  dplyr::summarise(                                        
    GO_terms = list(unique(GO_term)),                      
    .groups   = "drop"
  )

# gene → clusters lookup
gene2cluster <- ELAC_GO_elim_vs_WT_table %>%
  tidyr::unnest(                                          
    cols = unique_genes
  ) %>%
  dplyr::rename(                                         
    gene = unique_genes
  ) %>%
  dplyr::group_by(                                        
    gene
  ) %>%
  dplyr::summarise(                                       
    clusters = list(unique(cluster)),                     
    .groups  = "drop"
  )

#Join onto  DEG table
DEG_with_GO_and_cluster <- DEG_in_scRNA %>%
  dplyr::left_join(                                      
    gene2GO,
    by = "gene"
  ) %>%
  dplyr::left_join(                                      
    gene2cluster,
    by = "gene"
  ) %>%
  tidyr::replace_na(                                     
    list(
      GO_terms = list(character()),                      
      clusters = list(character())
    )
  )


#-------------------- Group the data and concatenate GO terms-------------------
grouped_GO <- ELAC_GO_elim_vs_WT %>%
  dplyr::group_by(cluster, dpa, Generalized.cell.population.group, expression) %>%
  dplyr::summarise(
    GO_terms = paste(GO_term, collapse = "; ")
  ) %>%
  dplyr::ungroup()

# Iterate over each row and print the formatted text
for(i in 1:nrow(grouped_GO)) {
  cat(
    "Cluster: ", grouped_GO$cluster[i], "\n",
    "Cell typ: ", grouped_GO$Generalized.cell.population.group[i], "\n",
    "Expresion: ", grouped_GO$expression[i], "\n",
    "DPA: ", grouped_GO$dpa[i], "\n",
    "GO terms: ", grouped_GO$GO_terms[i], "\n",
    "***\n",
    sep = ""
  )
}
#-------------------------------------------------------------------------------

head(ELAC_GO_elim_vs_WT)

MetabolicProcesses_data=ELAC_GO_elim_vs_WT[ELAC_GO_elim_vs_WT$cluster %in% c("Metabolic Processes"),] #Metabolism
MetabolicProcesses=unique(MetabolicProcesses_data$GO_term)
write.xlsx(ELAC_GO_elim_vs_WT,file="D:/Illumina_silencing/final_results/ELAC_GO_elim_vs_WT_clustered.xlsx")

BP_data1=MetabolicProcesses_data
#BP_data2=MetabolicProcesses_data[MetabolicProcesses_data$GO_term %in% MetabolicProcesses[c(91:length(MetabolicProcesses))],]

BP_data3=ELAC_GO_elim_vs_WT[ELAC_GO_elim_vs_WT$cluster %in% c("Developmental Processes","Cell Proliferation & Differentiation",
                                                              "Apoptosis & Cell Death","Specialized Structures"),] 
BP_data4=ELAC_GO_elim_vs_WT[ELAC_GO_elim_vs_WT$cluster %in% c("Immune Response","Response to Stimuli & Stress","Signaling Pathways"),] 
BP_data5=ELAC_GO_elim_vs_WT[ELAC_GO_elim_vs_WT$cluster %in% c("Binding Activities","Iron, Metal & Oxygen Binding ",
                                                              "Membrane Structures & Transport",
                                                              "Lipid Processes"),]
BP_data6=ELAC_GO_elim_vs_WT[ELAC_GO_elim_vs_WT$cluster %in% c("Organellar Regions & Activities",
                                                              "Protein Folding & Processing","Miscellaneous Activities"),]
BP_data7=ELAC_GO_elim_vs_WT[ELAC_GO_elim_vs_WT$cluster %in% c("Blood Related Processes","Neuronal Processes","Extracellular Regions"),]


ELAC_GO_elim_vs_WT$GO_term[ELAC_GO_elim_vs_WT$cluster=="Extracellular & Organellar Regions"]
ELAC_GO_elim_vs_WT$GO_term[ELAC_GO_elim_vs_WT$cluster=="Organellar Activities"]

BP_data1$cluster <- str_wrap(BP_data1$cluster, width = 10)
#BP_data2$cluster <- str_wrap(BP_data2$cluster, width = 10)
BP_data3$cluster <- str_wrap(BP_data3$cluster, width = 10)
BP_data4$cluster <- str_wrap(BP_data4$cluster, width = 10)
BP_data5$cluster <- str_wrap(BP_data5$cluster, width = 10)
BP_data6$cluster <- str_wrap(BP_data6$cluster, width = 10)
BP_data7$cluster <- str_wrap(BP_data7$cluster, width = 10)

length(unique(BP_data1$GO_term)) #41
length(unique(BP_data3$GO_term)) #27
length(unique(BP_data4$GO_term)) #21
length(unique(BP_data5$GO_term)) #25
length(unique(BP_data6$GO_term)) #18
length(unique(BP_data7$GO_term)) #26
rbind(BP_data1,BP_data3)
#

length(unique(BP_data5$GO_term)) #31
length(unique(BP_data6$GO_term)) #32
length(unique(BP_data7$GO_term)) #22
# unique(BP_data1$Significant)
# unique(BP_data2$Significant)
# unique(BP_data3$Significant)
# unique(BP_data4$Significant)
#unique(BP_data5$Significant)

all_BP_data=list(BP_data1,BP_data2,BP_data3,BP_data4,BP_data5,BP_data6,BP_data7)

colnames(BP_data1)

# Re-assigning the factors to ensure GO_term is reversed in order
for (i in 1:length(all_BP_data)) {
  all_BP_data[[i]]$GO_term <- factor(all_BP_data[[i]]$GO_term, levels = rev(sort(unique(all_BP_data[[i]]$GO_term))))
}
#setwd("D:/Illumina_silencing/final_results/GO_plots/new")

# Loop through each data frame in the list
for (i in seq_along(all_BP_data)) {
  df <- all_BP_data[[i]]
  cat("Unique values for DataFrame", i, ":\n")
  mycol=c("cluster","Generalized.cell.population.group")
  # Loop through each column in the data frame
  for (col_name in mycol) {
    unique_values <- unique(df[[col_name]])
    cat("Column:", col_name, "\n")
    print(unique_values)
    cat("\n")
  }
}



library(ggplot2)
library(dplyr)
library(stringr)

library(ggh4x)
#ls("package:ggh4x")
# Load required library
library(ggplot2)
library(ggh4x)




generate_temp_plot <- function(temp_df) {
  library(ggplot2)
  library(ggh4x)
  library(dplyr)
  
  # Group the data by 'cluster' and get the unique values of 'Generalized.cell.population.group'
  result <- temp_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      unique_groups = list(unique(Generalized.cell.population.group)),
      count_unique = n_distinct(Generalized.cell.population.group)
    ) %>%
    dplyr::ungroup()
  
  # Calculate generalized_levels and cluster_levels
  generalized_levels <- result %>%
    dplyr::summarise(total_unique_groups = sum(count_unique)) %>%
    pull(total_unique_groups)
  
  cluster_levels <- length(unique(temp_df$cluster))
  
  # Total number of strips on the y-axis
  total_y_strips <- generalized_levels + cluster_levels
  
  # Create the background_y list
  background_y_list <- c(
    rep(list(element_rect(fill = "grey85")), generalized_levels),  # For Generalized.cell.population.group levels
    rep(list(element_rect(fill = "grey70")), cluster_levels)       # For cluster levels
  )
  
  # Create the text_y list
  text_y_list <- c(
    rep(list(element_text(face = "bold", size = 10, angle = 0, margin = margin(r = 15, t = 35, b = 35))), generalized_levels),  # For Generalized.cell.population.group levels
    rep(list(element_text(face = "bold", size = 12, angle = 0)), cluster_levels)       # For cluster levels
  )
  
  # Number of cluster and generalized levels
  num_cluster_levels <- cluster_levels
  num_generalized_levels <- generalized_levels
  
  # Define the strip widths
  strip_width_values <- c(
    rep(unit(2, "cm"), num_cluster_levels),  # Width for 'cluster' strips
    rep(unit(0.5, "cm"), num_generalized_levels)  # Width for 'Generalized.cell.population.group' strips
  )
  
  # Create the background_y list (with outer and inner strips)
  background_y_list <- c(
    rep(list(element_rect(fill = "grey70")), num_cluster_levels),         # For cluster levels (outer strips)
    rep(list(element_rect(fill = "grey85")), num_generalized_levels)      # For Generalized.cell.population.group levels (inner strips)
  )
  
  # Adjust the text_y_list 
  text_y_list <- c(
    rep(list(element_text(face = "bold", size = 10, angle = -90, margin = margin(t = 5, b = 5, l = 2, r = 2))), num_cluster_levels),
    rep(list(element_text(face = "bold", size = 6, angle = 0, margin = margin(t = 5, b = 5, l = 2, r = 2))), num_generalized_levels)
  )
  
  # Generate the plot
  temp_plot <- ggplot(temp_df) +
    geom_point(
      aes(x = dpa, y = GO_term, color = Fisher.elim, size = Significant),
      alpha = 0.5
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 10),
      axis.text.x = element_text(
        angle = 0, vjust = 1, hjust = 0.5, size = 14
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(hjust = 0.5),
      panel.spacing = unit(0.1, "lines"),
      strip.placement = "outside"
    ) +
    scale_x_discrete(expand = c(0.01, 1)) +
    facet_nested(
      cluster + Generalized.cell.population.group ~ expression,
      scales = "free",
      space = "free",
      strip = strip_nested(
        # Set different background 
        background_x = list(
          element_rect(fill = "darkred"),
          element_rect(fill = "blue")
        ), 
        background_y = background_y_list,

        text_x = element_text(
          angle = 0, face = "bold", size = 12, margin = margin(t = 5, b = 5),
          color = "white"
        ),
        text_y = text_y_list
      )
    ) +
    labs(
      title = "Gene ontology for ELAC KD worms vs Control",
      x = "Comparison",
      y = "",
      color = "p-value",
      size = "Number of expressed genes"
    ) +
    scale_color_gradient2(
      low = "darkgreen",
      mid = "blue",
      high = "red"
    ) +
    scale_size_continuous(
      range = c(2, 6),
      breaks = c(2, 4, 6)
    )

  return(temp_plot)
}

head(BP_data1)
Metabolic_processess=generate_temp_plot(BP_data1)
Other_processess1=generate_temp_plot(rbind(BP_data3,BP_data4))
Other_processess2=generate_temp_plot(rbind(BP_data6,BP_data7))
Other_processess3=generate_temp_plot(BP_data5)
# Display the plot
print(Metabolic_processess) #Metabolic_processess
print(Other_processess1) #Other_processess1
print(Other_processess2) #Other_processess2

#setwd("D:/Illumina_silencing/final_results/GO_plots/with_cell_type_new")
# Save the plot with increased dimensions
ggsave("Metabolic_processess.pdf", Metabolic_processess, width = 10,height = 10)
ggsave("Other_processess1.pdf", Other_processess1, width = 10,height = 13)
ggsave("Other_processess2.pdf", Other_processess2, width = 10,height = 16)
ggsave("Other_processess3.pdf", Other_processess3, width = 10,height = 11)



library(reshape2)
temp_df=ELAC_GO_elim_vs_WT
library(data.table)
colnames(temp_df)
temp_df_long <- melt(setDT(temp_df[,c(3,8,11,12,14,17:32)]), id.vars = c("Significant","expression","Expressed_genes",
                                                                         "GO_term","dpa","cluster","Generalized.cell.population.group"), variable.name = "cell_population")


temp_df_long=temp_df_long[,-8]
head(temp_df_long)
colnames(temp_df_long)
library(dplyr)
library(tidyr)

# Split the Expressed_genes column by space (assuming space separates the gene IDs)
# You can also use another delimiter if necessary, e.g., comma (",")
temp_df_long_new <- temp_df_long %>%
  separate_rows(Expressed_genes, sep = " ") %>%  # Split the genes into individual rows
  dplyr::group_by(expression, dpa, cluster, cell_population) %>%
  dplyr::summarise(
    unique_genes = n_distinct(Expressed_genes),  # Count the unique genes
    Expressed_genes = paste(unique(Expressed_genes), collapse = " ")  # Concatenate unique genes into a single string
  ) %>%
  dplyr::ungroup()

# View the resulting summarized data
print(head(temp_df_long_new))



colnames(temp_df)
heatmap_data <- dcast(temp_df, cluster ~ Generalized.cell.population.group, value.var = "expression")
heatmap_data_matrix <- as.matrix(heatmap_data[,-1])  # Remove the GO_term column for matrix
rownames(heatmap_data_matrix) <- heatmap_data$cluster

heatmap(heatmap_data_matrix, scale = "row", col = colorRampPalette(c("blue", "white", "red"))(100), 
        main = "Gene Expression across Cell Populations")



ggplot(temp_df %>% 
         group_by(cluster) %>% 
         tally(),  # Use tally() to count the number of rows per group
       aes(x = reorder(cluster, -n), y = n)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(
    title = "Number of Significant GO Terms per Cluster",
    x = "Cluster", 
    y = "Number of GO Terms"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



library(igraph)

# Select the relevant columns and remove duplicate edges
head(temp_df[,c(3,8,11,12,14,17:32)])
library(dplyr)
library(igraph)

# Get the unique dpa values
dpa_values <- unique(temp_df$dpa)

# Loop over each dpa value and generate a network plot
for (dpa_value in dpa_values) {
  # Filter the data for the current dpa
  edges_dpa <- temp_df %>%
    filter(dpa == dpa_value) %>%
    dplyr::select(cluster, Generalized.cell.population.group) %>%
    distinct()  # Remove duplicate edges
  
  # Create a graph object for this dpa
  g_dpa <- graph_from_data_frame(edges_dpa)
  
  # Plot the network for this dpa
  plot(g_dpa, 
       vertex.label = V(g_dpa)$name, 
       vertex.color = "lightblue", 
       vertex.size = 5, 
       edge.arrow.size = 0.5, 
       main = paste("Network of GO clusters and cell populations for", dpa_value))
}


ggplot(temp_df, aes(x = cluster, y = Fisher.elim)) +
  geom_violin(aes(fill = cluster)) +
  scale_y_log10() +  # Log-transform the y-axis to better visualize the distribution
  labs(
    title = "Distribution of Fisher's p-values across Clusters",
    x = "Cluster", 
    y = "Fisher's p-value (log scale)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(temp_df, aes(x = Generalized.cell.population.group, fill = cluster)) +
  geom_bar(position = "stack") +
  labs(
    title = "GO Terms by Cell Population and Cluster",
    x = "Cell Population", 
    y = "Number of GO Terms"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

