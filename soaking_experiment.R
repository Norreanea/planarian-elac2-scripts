library(rstatix)
library(datarium)
library(installr)
#remotes::install_github("GegznaV/biostat")
library(biostat)
library(nlme)
library(boot)
#library(reCondition)
library(FSA)
library(chisq.posthoc.test)
library(RVAideMemoire)
library(tidyr)
#install.packages("psych")
library(ggstatsplot)
library(psych)
library(DescTools)
library(ggplot2)
library(dplyr)
library(ggpubr)
# Load the packages
library(multcomp)
library(broom)

pacman::p_load(readxl,dplyr,nptest,tidyverse,ggpubr,rstatix,datarium,car,ggplot2,ggpubr,lme4,
               lmerTest,emmeans,Matrix)
#?readxl::read_xlsx
################################################################################
################################################################################
#tRF level

# # Factors:
# # 
# #   Gene: 3 levels (ELAC2 KD, GFP Mock, WT)
# # Solution: 2 levels (planarian water, 5'tiRNAGlyGCC)
# # Groups: 6 combinations (e.g., ELAC2 KD + planarian water, GFP Mock + 5'tiRNAGlyGCC, etc.)
# # 
# # Replicates: 3 independent replicates per group.


soaking_experiment <- readxl::read_xlsx("D:/Elac2/final_results/tables/soaking_experiment_tRF_level.xlsx",
                                        sheet = "Sheet1")
soaking_experiment$Expression <- as.numeric(soaking_experiment$Expression)
unique(soaking_experiment$Condition)
# Perform Shapiro-Wilk test for each group
#shapiro.test(soaking_experiment$Percentage[soaking_experiment$Condition=="WT"])
(shapiro_results <- soaking_experiment[!(soaking_experiment$Condition %in% c("WT planarian water","GFP Mock no soak",
                                                                             "WT no soak","Elac2 KD planarian water")) ,] %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarize(
      shapiro_p = shapiro.test(Expression)$p.value,
      .groups = 'drop'
    )) 


summary_soaking_experiment <- soaking_experiment %>%
  dplyr::group_by(Condition,Group) %>% 
  dplyr::summarise(
    mean_exp = median(Expression),
    sd_exp = sd(Expression),
    n = n(),
    se_exp = sd_exp / sqrt(n),
    .groups = 'drop'
  ) %>%
  dplyr::mutate(
    Condition_num = as.numeric(factor(Condition))
  )

# Reorder 'Condition' levels so WT Water is first
summary_soaking_experiment$Condition <- factor(
  summary_soaking_experiment$Condition, 
  levels = c("WT planarian water", 
             unique(summary_soaking_experiment$Condition[summary_soaking_experiment$Condition != "WT planarian water"]))
)

# Perform Tukey HSD and filter for WT Water comparisons
?tukey_hsd
stat.test <- soaking_experiment %>%
  dplyr::group_by(Group) %>%
  tukey_hsd(Expression ~ Condition) %>%
  adjust_pvalue() %>%
  add_significance("p.adj") 

# stat.testWT <- stat.test %>%
#   dplyr::filter(group2 == "WT Water")
unique(stat.test$group1)
unique(stat.test$group2)

stat.test$group1 <- ifelse(stat.test$group1=="Elac2 KD 5'tiRNAGlyGCC","Elac2 KD 5'tiRNA Gly-GCC",stat.test$group1)
stat.test$group2 <- ifelse(stat.test$group2=="GFP Mock 5'tiRNAGlyGCC","GFP Mock 5'tiRNA Gly-GCC",stat.test$group2)
stat.test$group1 <- ifelse(stat.test$group1=="GFP Mock 5'tiRNAGlyGCC","GFP Mock 5'tiRNA Gly-GCC",stat.test$group1)
stat.test$group1 <- ifelse(stat.test$group1=="WT 5'tiRNAGlyGCC","WT 5'tiRNA Gly-GCC",stat.test$group1)


summary_soaking_experiment$Condition <- as.character(summary_soaking_experiment$Condition)
summary_soaking_experiment$Condition <- ifelse(summary_soaking_experiment$Condition=="Elac2 KD 5'tiRNAGlyGCC",
                                               "Elac2 KD 5'tiRNA Gly-GCC",
                                               summary_soaking_experiment$Condition)
summary_soaking_experiment$Condition <- ifelse(summary_soaking_experiment$Condition=="GFP Mock 5'tiRNAGlyGCC",
                                               "GFP Mock 5'tiRNA Gly-GCC",
                                               summary_soaking_experiment$Condition)
summary_soaking_experiment$Condition <- ifelse(summary_soaking_experiment$Condition=="WT 5'tiRNAGlyGCC",
                                               "WT 5'tiRNA Gly-GCC",
                                               summary_soaking_experiment$Condition)
# Create a comparison_id for each row in stat.testWT
stat.test_selected <- stat.test %>%
  dplyr::mutate(comparison_id = paste(Group, group1, group2, sep = "_vs_"))

comparison_data <- stat.test_selected %>%
  dplyr::rowwise() %>%
  dplyr::do({
    Group <- .$Group
    g1 <- .$group1
    g2 <- .$group2
    cid <- .$comparison_id
    
    d1 <- summary_soaking_experiment %>%
      dplyr::filter(Group == Group, Condition == g1)
    d2 <- summary_soaking_experiment %>%
      dplyr::filter(Group == Group, Condition == g2)
    # Bind test Condition and WT Water
    dplyr::bind_rows(d1, d2) %>%
      dplyr::mutate(comparison_id = cid)
  }) %>%
  dplyr::ungroup()


annot_data <- stat.test_selected %>%
  dplyr::mutate(y.position = 1.1 * max(comparison_data$mean_exp)) 


library(RColorBrewer)
display.brewer.pal(n = 10, name = "Set3")
# Retrieve the Set3 colors
set3_colors <- brewer.pal(n = 10, name = "Set3")
# Display the retrieved colors
print(set3_colors)

unique(comparison_data$Condition)[1]

wrapped_labels <- c(
  "Elac2 KD \n5'tiRNA Gly-GCC",
  "Elac2 KD \nplanarian water",
  "Elac2 KD \nscrambled",
  "GFP Mock \n5'tiRNA Gly-GCC",
  "GFP Mock \nplanarian water",
  "GFP Mock \nscrambled",
  "WT 5'tiRNA Gly-GCC",
  "WT planarian water",
  "WT scrambled" 
)
unique(comparison_data$Condition)
comparison_data<-na.omit(comparison_data)
comparison_data$Gene<-ifelse(startsWith(comparison_data$Condition,"WT"),"WT","GFP Mock")
comparison_data$Gene<-ifelse(startsWith(comparison_data$Condition,"ELAC2 KD"),"ELAC2 KD",comparison_data$Gene)



annot_data$y.position[10]=1.74
annot_data$y.position[11]=1.84
annot_data$y.position[12]=1.94
annot_data$y.position[13]=1.34
annot_data$y.position[14]=1.44
annot_data$y.position[15]=1.54

annot_data$y.position[1]=2.74
annot_data$y.position[2]=2.84
annot_data$y.position[3]=2.64
annot_data$y.position[4]=1.34
annot_data$y.position[5]=1.44
annot_data$y.position[6]=1.54
annot_data$y.position[7]=1.94
annot_data$y.position[8]=2.06
annot_data$y.position[9]=1.84

tRF_soaking_plot1 <- ggplot(
  soaking_experiment %>% filter(Group == "group4"),
  aes(x     = Condition,
      y     = Expression,
      fill  = Condition)   # or fill = Gene, whatever grouping you like
) +
  geom_boxplot(width = .6, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = .15),
              size     = 2,
              alpha    = .7) +
  stat_pvalue_manual(
    annot_data %>% filter(Group == "group4"),
    label      = "p.adj.signif",
    y.position = "y.position",
    tip.length = 0.01,
    inherit.aes= FALSE
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )+
  scale_x_discrete(labels = c("5` tiRNA \nGLY-GCC",
                              "5` tiRNA \nGLY-GCC",
                              "no RNA")) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, .05))) +
  labs(y = "Relative 5'-tiRNA Gly-GCC accumulation") 

tRF_soaking_plot1


tRF_soaking_plot2 <- ggplot(
  soaking_experiment %>% filter(Group == "group5"),
  aes(x     = Condition,
      y     = Expression,
      fill  = Condition)   # or fill = Gene, whatever grouping you like
) +
  geom_boxplot(width = .6, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = .15),
              size     = 2,
              alpha    = .7) +
  stat_pvalue_manual(
    annot_data %>% filter(Group == "group5"),
    label      = "p.adj.signif",
    y.position = "y.position",
    tip.length = 0.01,
    inherit.aes= FALSE
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )+
  scale_x_discrete(labels = c("no RNA",
                              "no RNA",
                              "no RNA")) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, .05))) +
  labs(y = "Relative 5'-tiRNA Gly-GCC accumulation") 

tRF_soaking_plot2

my_colors <- c(
  "ELAC2 KD 5'tiRNA Gly-GCC" = "#F8766D",
  "ELAC2 KD planarian water" = "#F8766D",
  "ELAC2 KD scrambled"       = "#F8766D",
  "GFP Mock 5'tiRNA Gly-GCC" = "#00BA38",
  "GFP Mock planarian water" = "#00BA38",
  "GFP Mock scrambled"       = "#00BA38",
  "WT 5'tiRNA Gly-GCC"       = "#619CFF",
  "WT planarian water"       = "#619CFF",
  "WT scrambled"             = "#619CFF"
)


soaking_experiment %>% filter(Group == "group4")
soaking_experiment %>% filter(Group == "group5")
library(lme4)  # For mixed-effects models
?lme4
# Fit the mixed-effects model
#unique(soaking_experiment$Condition)
soaking_experiment4 <- soaking_experiment %>% filter(Group == "group4")
soaking_experiment5 <- soaking_experiment %>% filter(Group == "group5")
model <- lmer(Expression ~ Condition + (1 | Replicate), data = soaking_experiment5)
# Perform Tukey's Honest Significant Difference (HSD) test
glht_result <- glht(model, linfct = mcp(Condition = "Tukey"))
# Summarize the results
summary_glht <- summary(glht_result)

#?glht
# Extract residuals and fitted values
residuals <- resid(model)
fitted <- fitted(model)

print(summary_glht)
# Tidy the glht result
tidy_glht <- broom::tidy(glht_result)
tidy_glht <- tidy_glht %>%
  separate(contrast, into = c("group1", "group2"), sep = " - ") 
# Add significance labels
(tidy_glht <- tidy_glht %>%
    add_significance("adj.p.value"))

########################################################################
soaking_experiment <- readxl::read_xlsx("D:/ELAC2/final_results/tables/soaking_experiment_gene_level.xlsx",
                                        sheet = "Sheet1")
unique(soaking_experiment$Condition)

as.data.frame(soaking_experiment)
soaking_experiment <- soaking_experiment[soaking_experiment$Target=="SMESG000048842.1",]

(shapiro_results <- soaking_experiment[soaking_experiment$Condition!="WT planarian water",] %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarize(
      shapiro_p = shapiro.test(Expression)$p.value,
      .groups = 'drop'
    )) 


# Perform Levene's Test
(levene_result <- leveneTest(Expression ~ Condition, data = soaking_experiment)) # The variances are considered equal across groups

soaking_experiment
unique(soaking_experiment$Condition)

summary_soaking_experiment <- soaking_experiment %>%
  dplyr::group_by(Condition,Target) %>% 
  dplyr::summarise(
    mean_exp = median(Expression),
    sd_exp = sd(Expression),
    n = n(),
    se_exp = sd_exp / sqrt(n),
    .groups = 'drop'
  ) %>%
  dplyr::mutate(
    Condition_num = as.numeric(factor(Condition))
  )


# Reorder 'Condition' levels so WT Water is first
summary_soaking_experiment$Condition <- factor(
  summary_soaking_experiment$Condition, 
  levels = c("WT planarian water", unique(summary_soaking_experiment$Condition[summary_soaking_experiment$Condition != "WT planarian water"]))
)


soaking_experiment$Condition<-as.factor(soaking_experiment$Condition)

# Fit the linear model
library(lme4)  # For mixed-effects models
?lme4
# Fit the mixed-effects model
unique(soaking_experiment$Condition)
model <- lmer(Expression ~ Condition + (1 | Replicate), data = soaking_experiment)
# Perform Tukey's Honest Significant Difference (HSD) test
glht_result <- glht(model, linfct = mcp(Condition = "Tukey"))
# Summarize the results
summary_glht <- summary(glht_result)

?glht
# Extract residuals and fitted values
residuals <- resid(model)
fitted <- fitted(model)

# Plot
ggplot(data = NULL, aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs. Fitted Values",
       x = "Fitted Values",
       y = "Residuals") +
  theme_minimal()
shapiro.test(residuals)



print(summary_glht)
# Tidy the glht result
tidy_glht <- broom::tidy(glht_result)
tidy_glht <- tidy_glht %>%
  separate(contrast, into = c("group1", "group2"), sep = " - ") 
# Add significance labels
(tidy_glht <- tidy_glht %>%
    add_significance("adj.p.value"))

stat.test<-tidy_glht
unique(stat.test$group1)
unique(stat.test$group2)

stat.test$group1 <- ifelse(stat.test$group1=="ELAC2 KD 5'tiRNAGlyGCC","ELAC2 KD 5'tiRNA Gly-GCC",stat.test$group1)
stat.test$group2 <- ifelse(stat.test$group2=="ELAC2 KD 5'tiRNAGlyGCC","ELAC2 KD 5'tiRNA Gly-GCC",stat.test$group2)
stat.test$group2 <- ifelse(stat.test$group2=="GFP Mock 5'tiRNAGlyGCC","GFP Mock 5'tiRNA Gly-GCC",stat.test$group2)
stat.test$group1 <- ifelse(stat.test$group1=="GFP Mock 5'tiRNAGlyGCC","GFP Mock 5'tiRNA Gly-GCC",stat.test$group1)
stat.test$group1 <- ifelse(stat.test$group1=="WT 5'tiRNAGlyGCC","WT 5'tiRNA Gly-GCC",stat.test$group1)
stat.test$group2 <- ifelse(stat.test$group2=="WT 5'tiRNAGlyGCC","WT 5'tiRNA Gly-GCC",stat.test$group2)

stat.test_selected <- stat.test
as.data.frame(stat.test_selected)
as.data.frame(stat.test)

summary_soaking_experiment$Condition <- as.character(summary_soaking_experiment$Condition)
summary_soaking_experiment$Condition <- ifelse(summary_soaking_experiment$Condition=="ELAC2 KD 5'tiRNAGlyGCC","ELAC2 KD 5'tiRNA Gly-GCC",
                                               summary_soaking_experiment$Condition)
summary_soaking_experiment$Condition <- ifelse(summary_soaking_experiment$Condition=="GFP Mock 5'tiRNAGlyGCC","GFP Mock 5'tiRNA Gly-GCC",
                                               summary_soaking_experiment$Condition)
summary_soaking_experiment$Condition <- ifelse(summary_soaking_experiment$Condition=="WT 5'tiRNAGlyGCC","WT 5'tiRNA Gly-GCC",
                                               summary_soaking_experiment$Condition)

# # Create a comparison_id for each row in stat.testWT
stat.test_selected <- stat.test_selected %>%
  dplyr::mutate(comparison_id = paste(group1, group2, sep = "_vs_"))


# Plot each comparison separately

annot_data <- stat.test_selected %>%
  dplyr::mutate(y.position = 1.1 * max(summary_soaking_experiment$mean_exp)) 



# Define a custom labeller function to extract the gene identifier
gene_labeller <- function(comparison_id) {
  # Split the comparison_id by "_vs_" and take the first part
  gene_id <- sapply(str_split(comparison_id, "_vs_"), `[`, 1)
  return(gene_id)
}

gene_labeller <- labeller(comparison_id = function(x) str_extract(x, "^[^_]+"))


my_colors <- c("ELAC2 KD" = "#F8766D", "GFP Mock" = "#00BA38", "WT" = "#619CFF")
unique(comparison_data$comparison_id)
my_colors = c(
  "ELAC2 KD 5'tiRNAGlyGCC" = "#F8766D",
  "ELAC2 KD planarian water" = "#F8766D",
  "ELAC2 KD scrambled" = "#F8766D",
  "GFP Mock 5'tiRNAGlyGCC" = "#00BA38",
  "GFP Mock planarian water" = "#00BA38",
  "GFP Mock scrambled" = "#00BA38",
  "WT 5'tiRNAGlyGCC" = "#619CFF",
  "WT planarian water"= "#619CFF",
  "WT scrambled" = "#619CFF"
  
)
my_colors_df=as.data.frame(my_colors)
colnames(my_colors_df)="color"
my_colors_df$Condition = rownames(my_colors_df)
comparison_data = merge(summary_soaking_experiment,my_colors_df, by = "Condition", all = TRUE)
comparison_data = na.omit(comparison_data)
as.data.frame(annot_data)
which(annot_data$adj.p.value.signif!="ns")
annot_data$y.position[1] <- 2.32
annot_data$y.position[2] <- 2.42
annot_data$y.position[10] <- 2.52
annot_data$y.position[12] <- 2.62
annot_data$y.position[14] <- 2.72
annot_data$y.position[15] <- 2.82
annot_data$y.position[16] <- 2.92
annot_data$y.position[17] <- 3.02
annot_data$y.position[18] <- 3.12
annot_data$y.position[20] <- 3.22
annot_data$y.position[21] <- 3.32
#SMESG000048842.1
#as.data.frame(annot_data[annot_data$Target=="SMESG000048842.1",])
unique(comparison_data$Condition)
wrapped_labels <- c(
  "ELAC2 KD \n5'tiRNA Gly-GCC",
  "ELAC2 KD \nplanarian water",
  "ELAC2 KD \nscrambled",
  "GFP Mock \n5'tiRNA Gly-GCC",
  "GFP Mock \nplanarian water",
  "GFP Mock \nscrambled",
  "WT 5'tiRNA Gly-GCC",
  "WT planarian water",
  "WT scrambled" 
)

# annot_data$y.position[2] <- 1.52
# annot_data$y.position[3] <- 1.72
as.data.frame(annot_data)
SMESG000048842_plot=ggplot(comparison_data[comparison_data$Target=="SMESG000048842.1",], aes(x = Condition, y = mean_exp, fill = Condition)) +
  geom_col(position = position_dodge(width=0.8)) +
  geom_errorbar(aes(ymin = mean_exp - se_exp, ymax = mean_exp + se_exp),
                width = 0.2, position = position_dodge(width=0.8)) +
  #facet_wrap(~ comparison_id, scales = "free_x", ncol = 6, labeller = gene_labeller,as.table = FALSE) +
  stat_pvalue_manual(
    #annot_data,
    annot_data[annot_data$adj.p.value.signif!="ns",],
    #annot_data[annot_data$p.adj.signif!="ns" & annot_data$comparison_id %in% selection,],
    label = "adj.p.value.signif",
    y.position = "y.position", # the position for annotation
    tip.length = 0.01,
    inherit.aes = FALSE
  ) +
  labs(
    title = "SMESG000048842.1",
    y = "Relative gene expression"
  ) +
  scale_fill_manual(values = my_colors) +
  theme_classic(base_size = 12) +
  #guides(fill="none") +
  theme(#axis.text.x = element_text(angle=45, hjust=1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 8))

as.data.frame(annot_data)
comparison_data$Gene<-ifelse(startsWith(comparison_data$Condition,"WT"),"WT","GFP Mock")
comparison_data$Gene<-ifelse(startsWith(comparison_data$Condition,"ELAC2 KD"),"ELAC2 KD",comparison_data$Gene)
selected_annot_data=annot_data
selected_annot_data=selected_annot_data[selected_annot_data$group1=="ELAC2 KD 5'tiRNA Gly-GCC" |
                                          selected_annot_data$group2=="ELAC2 KD 5'tiRNA Gly-GCC",]
selected_annot_data=selected_annot_data[selected_annot_data$adj.p.value.signif!="ns" |
                                          selected_annot_data$group1=="WT planarian water",]
selected_annot_data$y.position[3]=2.52
(SMESG000048842_plot_another=ggplot(comparison_data[comparison_data$Target=="SMESG000048842.1",], 
                                    aes(x = Condition, y = mean_exp, fill = Gene)) +
    geom_col(position = position_dodge(width=0.8)) +
    geom_errorbar(aes(ymin = mean_exp - se_exp, ymax = mean_exp + se_exp),
                  width = 0.2, position = position_dodge(width=0.8)) +
    #facet_wrap(~ comparison_id, scales = "free_x", ncol = 6, labeller = gene_labeller,as.table = FALSE) +
    stat_pvalue_manual(
      #annot_data,
      selected_annot_data,
      #annot_data[annot_data$adj.p.value.signif!="ns",],
      #annot_data[annot_data$p.adj.signif!="ns" & annot_data$comparison_id %in% selection,],
      label = "adj.p.value.signif",
      y.position = "y.position", # the position for annotation
      tip.length = 0.01,
      inherit.aes = FALSE
    ) +
    labs(
      title = "SMESG000048842.1",
      y = "Relative SMESG000048842.1 expression"
    ) +
    scale_x_discrete(labels = c(rep(c("5` tRNA half \nGLY-GCC","no RNA","scrambled \nRNA"),3))) +
    #scale_fill_manual(values = my_colors) +
    theme_classic(base_size = 12) +
    #guides(fill="none") +
    theme(#axis.text.x = element_text(angle=90, hjust=1),
      #axis.text.x = element_blank(),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5)))
# legend.title = element_text(size = 8))


library(dplyr)
library(ggplot2)
library(ggpubr)      # for stat_pvalue_manual()

# Pull out the raw replicates for your gene of interest
df_gene <- soaking_experiment %>%
  dplyr::filter(Target == "SMESG000048842.1")

#  Make sure 'Condition' is a factor in the order you want on the x-axis
df_gene$Condition <- factor(df_gene$Condition,
                            levels = c(
                              "ELAC2 KD 5'tiRNAGlyGCC",
                              "ELAC2 KD planarian water",
                              "ELAC2 KD scrambled",
                              "GFP Mock 5'tiRNAGlyGCC",
                              "GFP Mock planarian water",
                              "GFP Mock scrambled",
                              "WT 5'tiRNAGlyGCC",
                              "WT planarian water",
                              "WT scrambled"
                            )
)

# Subset annotation table to only those comparisons for this gene
ann <- annot_data %>%
  filter(
    TRUE, 
    adj.p.value.signif != "ns"
  ) %>%
  # the `x` / `xend` mappings come from  group1/group2 columns
  mutate(
    xmin = as.numeric(factor(group1, levels=levels(df_gene$Condition))),
    xmax = as.numeric(factor(group2, levels=levels(df_gene$Condition)))
  )

# 4. Plot
tRF_soaking_plot3 <- ggplot(df_gene, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot(width = .6, outlier.shape = NA) +
  geom_jitter(width = .15, size = 2, alpha = .7) +
  stat_pvalue_manual(
    ann,
    label        = "adj.p.value.signif",
    xmin         = "xmin", 
    xmax         = "xmax",
    y.position   = "y.position",
    tip.length   = .01,
    inherit.aes  = FALSE
  ) +
  scale_fill_manual(values = my_colors) +
  #scale_fill_brewer(palette = "Set2") +   # or scale_fill_manual(values = my_colors)
  scale_x_discrete(labels = c(
    "5` tiRNA\nGLY-GCC", "no RNA", "scrambled\nRNA",
    "5` tiRNA\nGLY-GCC", "no RNA", "scrambled\nRNA",
    "5` tiRNA\nGLY-GCC", "no RNA", "scrambled\nRNA"
  )) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, .05))
  ) +
  labs(
    title = "",
    y     = "Relative SMESG000048842.1 expression",
    x     = ""
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

#tRF_soaking_plot1,tRF_soaking_plot2,tRF_soaking_plot3


ggarrange(ggarrange(tRF_soaking_plot1,tRF_soaking_plot2,ncol =2, nrow=1),
          tRF_soaking_plot3,nrow = 2,ncol=1)






unique(comparison_data$comparison_id)

p1 <- ggplot() + theme_void()
# ggarrange(ggarrange(p1,SMESG000048842_plot_another,nrow=1,legend="none",labels=c("A","B")),
#           tRF_soaking_plot,nrow=2,common.legend=TRUE,legend="bottom",labels = "C")
ggarrange(ggarrange(p1,tRF_soaking_plot1,tRF_soaking_plot2,nrow=1,legend="none",labels=c("A","B")),
          SMESG000048842_plot_another,nrow=2,common.legend=TRUE,legend="top",labels = "C")


ggarrange(tRF_soaking_plot1,tRF_soaking_plot1,ncol =2)


##Supplementary figure for tRNA half soaking
ggarrange(tRF_soaking_plot3,tRF_soaking_plot4,tRF_soaking_plot5,ncol =3,common.legend=FALSE,legend="top",labels=c("A","B","C"))

