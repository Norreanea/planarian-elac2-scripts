#Load libraries
if (!require("pacman")) install.packages("pacman") 
pacman::p_load(readxl,dplyr,nptest,tidyverse,ggpubr,rstatix,datarium,car,ggplot2,ggpubr)

#Read file
(all_chimeras=read_excel("G:/General/documents/AZ/chimeras_PCR.xlsx"))
# reg_chimeras=unique(all_chimeras$Region)
# type_chimeras=unique(all_chimeras$Type)
# gene_chimeras=unique(all_chimeras$Gene)
all_chimeras$gene_reg_chimeras=paste(all_chimeras$Gene,all_chimeras$Region,sep="_")
gene_reg_chimeras=unique(all_chimeras$gene_reg_chimeras)
all_chimeras$type_reg_chimeras=paste(all_chimeras$Type,all_chimeras$Region,sep="_")
type_reg_chimeras=unique(all_chimeras$type_reg_chimeras)

# data=all_chimeras
# print(data,n=100)
all_chimeras$group=paste(all_chimeras$Region,all_chimeras$Type,sep="_")
all_chimeras$group_gene=paste(all_chimeras$group,all_chimeras$Gene,sep="_")
group_gene=unique(all_chimeras$group_gene)


all_chimeras[all_chimeras$gene_reg_chimeras==gene_reg_chimeras[1],]
for (i in 1:length(gene_reg_chimeras)) {
  reg_n=all_chimeras[all_chimeras$gene_reg_chimeras==gene_reg_chimeras[i],]
  print(paste("Levene’s test (to check if the variances of a variables  are equal) for ",gene_reg_chimeras[i],":"))
  result = leveneTest(FC ~ Replicate, reg_n)
  print(result)
  print(paste("---------------------------------------------------------------------------------"))
}


for (i in 1:length(type_reg_chimeras)) {
  reg_n=all_chimeras[all_chimeras$type_reg_chimeras==type_reg_chimeras[i],]
  print(paste("One way ANOVA for ",type_reg_chimeras[i],":"))
  btwn <- aov(FC ~ Gene, data=reg_n)
  summary(btwn)
  result=TukeyHSD(btwn)
  #result = leveneTest(FC ~ Replicate, reg_n)
  print(result)
  print(paste("---------------------------------------------------------------------------------"))
}
unique(all_chimeras$group_gene)
all_chimeras_no_GFP=all_chimeras[all_chimeras$Gene!="GFP",]
all_chimeras$gene_reg_type=paste(all_chimeras$gene_reg_chimeras,all_chimeras$group)
all_chimeras %>%
  group_by(group) %>%
  games_howell_test(FC ~ Gene)


library(PMCMRplus)
for (i in 1:length(type_reg_chimeras)) {
  reg_n=all_chimeras[all_chimeras$type_reg_chimeras==type_reg_chimeras[i],]
  print(paste("One way ANOVA for ",type_reg_chimeras[i],":"))
  btwn <- aov(FC ~ Gene, data=reg_n)
  summary(btwn)
  result=dunnettT3Test(btwn)
  #result = leveneTest(FC ~ Replicate, reg_n)
  print(result)
  print(paste("---------------------------------------------------------------------------------"))
}
?glm


for (i in 1:length(type_reg_chimeras)) {
  reg_n=all_chimeras[all_chimeras$type_reg_chimeras==type_reg_chimeras[i],]
  print(paste("GLM for ",type_reg_chimeras[i],":"))
  btwn=glm(formula = FC ~ Gene, data=reg_n,family = quasipoisson())
  #btwn <- t.test(FC ~ Gene, data=reg_n,var.equal = FALSE)
  summary(btwn)
  #result=TukeyHSD(btwn)
  #result = leveneTest(FC ~ Replicate, reg_n)
  print(result)
  print(paste("---------------------------------------------------------------------------------"))
}

reg_n=all_chimeras_no_GFP[all_chimeras_no_GFP$type_reg_chimeras==type_reg_chimeras[1],]
print(paste("t-Student test for ",type_reg_chimeras[1],":"))
(btwn <- t.test(FC ~ Gene, data=reg_n,var.equal = TRUE,alternative ="greater"))
?t.test
?TukeyHSD
for (i in 1:length(type_reg_chimeras)) {
  reg_n=all_chimeras_no_GFP[all_chimeras_no_GFP$type_reg_chimeras==type_reg_chimeras[i],]
  print(paste("t-Student test for ",type_reg_chimeras[i],":"))
  btwn <- t.test(FC ~ Gene, data=reg_n,var.equal = FALSE,alternative ="greater")
  print(btwn)
  #result=TukeyHSD(btwn)
  #result = leveneTest(FC ~ Replicate, reg_n)
  #print(result)
  print(paste("---------------------------------------------------------------------------------"))
}

# Kruskal-Wallis test to compare FC between genes within each group (region + type)
(kruskal_test_results = all_chimeras %>%
  dplyr::group_by(group) %>%
  dplyr::summarize(p_value = kruskal.test(FC, Gene)$p.value))
?wilcox.test
print(kruskal_test_results)
# Assuming your data is already loaded as 'all_chimeras'
# Make sure you have loaded the necessary packages: dplyr, broom

# # Pairwise Wilcoxon rank-sum test (Mann-Whitney U test) to compare FC between genes within each group (region + type)
pairwise_wilcoxon_test_results <- all_chimeras %>%
  group_by(Region, Type) %>%
  do(test_result = tidy(pairwise.wilcox.test(.$FC, .$Gene, p.adjust.method = "none")))

# # Unnest the results
pairwise_wilcoxon_test_results <- pairwise_wilcoxon_test_results %>%
  tidyr::unnest(cols = test_result)

# Print the results
print(pairwise_wilcoxon_test_results)


# Perform two-way ANOVA
anova_results <- all_chimeras %>%
  dplyr::group_by(Region, Type) %>%
  aov(FC ~ Gene:Region:Type, data = .)
tukey_results <- TukeyHSD(anova_results)
# Summarize the ANOVA results
summary(anova_results)
reg_n=all_chimeras[all_chimeras$type_reg_chimeras==type_reg_chimeras[1],]
btwn <- aov(FC ~ Gene, data=reg_n)
summary(btwn)
TukeyHSD_out<-TukeyHSD(anova_results)
TukeyHSD_out<-as.data.frame(TukeyHSD_out)
TukeyHSD_out[TukeyHSD_out$`Gene:Region:Type`$p.adj<0.05,]

# Filter groups with more than one gene for the Wilcoxon rank-sum test
filtered_data <- all_chimeras %>%
  group_by(Region, Type) %>%
  filter(n_distinct(Gene) >= 2)

# Check which groups have more than one unique gene
multi_gene_groups <- filtered_data %>%
  group_by(Region, Type) %>%
  filter(n_distinct(Gene) >= 2) %>%
  distinct(Region, Type)

# Wilcoxon rank-sum test (Mann-Whitney U test) to compare FC between genes within each group (region + type)
wilcoxon_test_results <- multi_gene_groups %>%
  group_split() %>%
  map(~ wilcox.test(FC ~ Gene, data = filter(filtered_data, Region == .$Region[1] & Type == .$Type[1])) %>%
        broom::tidy() %>%
        select(Region, Type, p.value))


# Aggregate to calculate the variance for each group
group_names = paste(data$Region, data$Type, data$Gene, sep = "-")
variance_results = aggregate(FC ~ group_names, data = data, var)

# Rename the columns for clarity
colnames(variance_results) = c("Group", "Variance")

# Print the variance results
print(variance_results)

# ggplot(all_chimeras, aes(y=FC, x=Gene, color = factor(Gene)))+
#   geom_boxplot()+
#   theme( legend.position = "none" ) + theme_bw()+facet_grid(~Region+Type)

# ggplot(all_chimeras, aes(y=FC, x=Gene, color = factor(Gene)))+
#   geom_bar()+
#   theme( legend.position = "none" ) + theme_bw()+facet_grid(~Region+Type)

# Box plot facetted by "dose"
?ggboxplot
###############################################################################
#PLOT
all_chimeras$Condition=all_chimeras$Gene
all_chimeras$Condition=ifelse(all_chimeras$Condition=="GFP","GFP mock",all_chimeras$Condition)
all_chimeras$Condition=ifelse(all_chimeras$Condition=="ELAC2","ELAC2 KD",all_chimeras$Condition)
all_chimeras$Type=str_to_sentence(all_chimeras$Type)
all_chimeras$Region=str_to_sentence(all_chimeras$Region)

?tukey_hsd()

# Statistical test
(stat.test <- all_chimeras %>%
  group_by(Region,Type) %>%
  tukey_hsd(FC ~ Condition) %>%
  adjust_pvalue() %>%
  add_significance("p.adj"))
stat.test[stat.test$p.adj.signif!="ns",]
?stat_pvalue_manual

(p = ggboxplot(all_chimeras, x = "Condition", y = "FC",
               width=0.2,size=1,
               facet.by = c("Region","Type"), short.panel.labs = TRUE,bxp.errorbar = TRUE)
  + stat_pvalue_manual(stat.test, label = "p.adj.signif",hide.ns=TRUE,y.position = c(7,5,5,6))
  +theme_bw(base_size = 22))+
  geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0),
             pch=21, aes(fill=Condition), size=3,show.legend = T)+  
  theme(text = element_text(size = 20),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background=element_rect(fill="white"),legend.position="none")+
  labs(y = "Fold change")+scale_x_discrete(labels = c("ELAC2\nKD","GFP\nmock","WT"))

#boxplot with no statistics
(p = ggboxplot(all_chimeras, x = "Condition", y = "FC",
               width=0.2,size=1,
               facet.by = c("Region","Type"), short.panel.labs = TRUE,bxp.errorbar = TRUE)+theme_bw(base_size = 22))+
    geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0),
               pch=21, aes(fill=Condition), size=3,show.legend = T)+  
  theme(text = element_text(size = 20),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background=element_rect(fill="white"),legend.position="none")+
  labs(y = "Fold change")+scale_x_discrete(labels = c("ELAC2\nKD","GFP\nmock","WT"))
p+ stat_pvalue_manual(stat.test, label = "p.adj", y.position = 3)
as.data.frame(stat.test)

?stat_compare_means
my_comparisons = list(c("ELAC2", "GFP"), c("ELAC2", "WT"), c("GFP", "WT") )
ggplot(all_chimeras, aes(y=FC, x=Gene, color = factor(Gene)))+
  geom_boxplot()+
  theme( legend.position = "none" ) + theme_bw()+facet_grid(~Region+Type)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 10,aes(label = paste0("p = ", after_stat(p.format))))     # Add global p-value


?stat_compare_means
