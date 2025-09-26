
#LIBRARIES
#------------------------------libraries----------------------------------------
# install.packages("dplyr", type = "source")
# install.packages("xfun", type = "source")
# remove.packages("rlang")
# install.packages("rlang")
# 
# update.packages()
# # Restart R session
# restartSession()
# 
# # Reinstall rlang
# install.packages("rlang")
# 
# # Install DESeq2
# install.packages("DESeq2")
# 
# install.packages("rlang", dependencies = TRUE)
# devtools::install_github("dzhang32/ggtranscript")
#remotes::install_github("dzhang32/ggtranscript")
# # install via CRAN
# install.package("ggcoverage")
# 
# # install via Github
# # install.package("remotes")   #In case you have not installed it.
# remotes::install_github("showteeth/ggcoverage")
# if (!require("pacman")) install.packages("pacman") 
# #install required libraries (if needed) and load them
pacman::p_load(fdrtool,tximport,tximportData,edgeR,DESeq2,readr,tidyverse,RColorBrewer,pheatmap,DEGreport,ggplot2,ggrepel,"vsn",goeveg,ape,rtracklayer,ggtree,
               pheatmap,tidyverse,hrbrthemes,viridis,data.table,reshape,gplots,ggpubr,VennDiagram,limma,WGCNA,sm,network,venn,xlsx,openxlsx ,
               dplyr,topGO,ALL,AnnotationDbi,org.Hs.eg.db,pathview,gage,gageData,plyr,strex,UpSetR,GOfuncR,gridExtra,grid,data.table,ggh4x,PCAtools)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "tidyverse", "pheatmap", "DEGreport", "ggplot2", "ggrepel", "vsn", "goeveg", "ggtree",
#                      "pheatmap", "tidyverse", "hrbrthemes", "viridis", "ggpubr", "WGCNA", "UpSetR", "gridExtra", "ggh4x", "PCAtools"),force =
#                        TRUE)
# install.packages("DESeq2")
# library(rlang)
# library(DESeq2)
# library(dplyr)
#devtools::install_github("dzhang32/ggtranscript")
#remotes::install_github("showteeth/ggcoverage")
library(ggtranscript)
#install.packages("ggtranscript")
#BiocManager::install("plotgardener")
#install.packages("BSgenome")
library(ggcoverage)
# mtDNA_anno_gff3=read.gff("E:/genome/mtDNA_anno.gff3", na.strings = c(".", "?"), GFF3 = TRUE)
# mtDNA_anno_gff3$attributes
# mtDNA_anno_gtf <- import("E:/genome/mtDNA_anno.gtf")
# 
# attr_gff3=mtDNA_anno_gff3$attributes
# typeof(attr_gff3)
# attr_gff3[1]
# mtDNA_anno_gtf$transcript_id
# 
# # Extract elements with "Parent" and "product"
# matches <- grep("Parent=rna.*;product=.*", attr_gff3, value = TRUE)
# 
# # Extract "Parent" and "product" values using regular expressions
# parent <- gsub(".*Parent=([^;]+);.*", "\\1", matches)
# product <- gsub(".*product=([^;]+)", "\\1", matches)
# 
# 
# # Create data frame with desired columns
# (features <- data.frame(Parent = parent, Product = product))
# 
# 
# features$Product[features$Parent==mtDNA_anno_gtf$transcript_id[3]]
# #change feature name in gtf
# for (i in 1:length(mtDNA_anno_gtf)){
#   mtDNA_anno_gtf$transcript_id[i]=ifelse(mtDNA_anno_gtf$transcript_id[i] %in% features$Parent, features$Product[features$Parent==mtDNA_anno_gtf$transcript_id[i]],mtDNA_anno_gtf$transcript_id[i])
# }
# 
# mtDNA_anno_gtf$gene_id=mtDNA_anno_gtf$transcript_id
# mtDNA_anno_gtf$gene_name=ifelse(is.na(mtDNA_anno_gtf$gene_name),mtDNA_anno_gtf$gene_id,mtDNA_anno_gtf$gene_name)
# export(mtDNA_anno_gtf, 'E:/genome/mtDNA_anno_corrected.gtf')
mtDNA_anno_gtf=rtracklayer::import('D:/genome/genome/mtDNA_anno_corrected.gtf')


#loading all genome and transcriptome
(cov_list=list.files("E:/Illumina/PARN_ELAC2_silencing/mRNA/mapped_to_chr_and_mtDNA/mtDNA_coverage/", pattern=NULL, all.files=FALSE,
                     full.names=FALSE))

new_name <- sub("_coverage\\.txt$", "", cov_list)
# 
# ELAC_L131=read.table("E:/Illumina/PARN_ELAC2_silencing/mRNA/mapped_to_chr_and_mtDNA/mtDNA_coverage/Elac_L13_batch1_coverage.txt",
#                      col.names = c("seqnames","start","end","score"))
# 
# ELAC_L131$Type="ELAC_L13"


# ?ggcoverage
# 
# (basic.coverage = ggcoverage(data = ELAC_L131,range.position = "out"))
# ggcoverage(data = ELAC_L131)
# 
# 
# 
# p17 <- autoplot(txdb, genesymbol["RBM17"])
# plotRangesLinkedToData(exon.new, stat.y = c("sample1", "sample2"), annotation = list(p17))
all_cov_df = data.frame()
for (j in 1:length(cov_list)) {
  print(paste("Proceeding file",j))
  sample_df=read.table(paste("E:/Illumina/PARN_ELAC2_silencing/mRNA/mapped_to_chr_and_mtDNA/mtDNA_coverage/",cov_list[j],
                             sep=""),col.names = c("seqnames","start","end","score"))
  sample_df$Type=new_name[j]
  all_cov_df=rbind(all_cov_df,sample_df)
  #assign(paste("trans_genome_RNA",j,sep=""),mapped_seq)
}

all_cov_df=all_cov_df[!startsWith(all_cov_df$Type,"PARN"),]
# 
# 
#     # Filter long read gtf/gff for gene of interest
#     # gtf_filtered <- gtf[gtf$gene_id == gene_id]
#     # 
#     # # loci used to filter data
#     # locus_subset <- GRanges(seqnames = seqnames, ranges = IRanges(start = start, end = end), strand = strand)
#     # 
#     # # coverage data
#     # coverage_data <- rtracklayer::import.bw(coverage) %>% subsetByOverlaps(., locus_subset) %>% as.data.frame()
#     
#     # Plot transcripts
#     exons <- data.frame(mtDNA_anno_gtf) %>% dplyr::filter(type == "exon")
#     introns <- exons %>% to_intron(group_var = "transcript_id")
#     CDS <- data.frame(mtDNA_anno_gtf) %>% dplyr::filter(type == "CDS")
#     
#     transcript_plot <-
#       exons %>%
#       ggplot(aes(
#         xstart = start, 
#         xend = end, 
#         y = transcript_id)) +
#       geom_range(fill = "white",
#                  height = 0.25) +
#       geom_range(data = CDS) +
#       geom_intron(
#         data = introns,
#         arrow.min.intron.length = 500,
#         arrow = grid::arrow(ends = "first", length = grid::unit(0.1, "inches"))
#       ) +
#       labs(y = "Transcript name", x = "") +
#       xlim(start(locus_subset), end(locus_subset))
#     
#     # coverage data
#     coverage_plot <-
#       coverage_data %>%
#       ggplot(aes(
#         xmin = start,
#         xmax = end,
#         ymin = 0,
#         ymax = score
#       )) +
#       geom_rect(show.legend = F, alpha = 0.8) +
#       xlim(start(locus_subset), end(locus_subset))
#     
#     # Final plot
#     transcript_coverage_plot <-
#       plot_grid(
#         transcript_plot,
#         coverage_plot,
#         ncol = 1,
#         align = "hv",
#         rel_heights = c(1, 3.5),
#         axis = "lr"
#       )
#     
#     return(transcript_coverage_plot)
# 
# 
# 
# TranscriptCoveragePlot(mtDNA_anno_gtf)

# 
# WT_L131=read.table("E:/Illumina/PARN_ELAC2_silencing/mRNA/mapped_to_chr_and_mtDNA/mtDNA_coverage/WT_L13_batch1_coverage.txt",
#                      col.names = c("seqnames","start","end","score"))
# 
# WT_L131$Type="WT_L13"
# 
# 
# 
# ELAC_L131
# ELAC_WT=rbind(ELAC_L131,WT_L131)

library(ggplot2)
#library(babynames) # provide the dataset: a dataframe called babynames
library(dplyr)
library(hrbrthemes)
library(viridis)



all_cov_df=all_cov_df[,-3]

ELAC_cov=subset(all_cov_df,startsWith(all_cov_df$Type,"Elac"))
WT_cov=subset(all_cov_df,startsWith(all_cov_df$Type,"WT"))

GFP_cov=subset(all_cov_df,startsWith(all_cov_df$Type,"GFP"))

library(tidyr)
?median

ELAC_cov_wide <- spread(ELAC_cov, Type, score)
ELAC_cov_wide$ELAC_L13=ELAC_cov_wide$Elac_L13_batch1+ELAC_cov_wide$Elac_L13_batch2
ELAC_cov_wide$ELAC_L23=ELAC_cov_wide$Elac_L23_batch1+ELAC_cov_wide$Elac_L23_batch2
ELAC_cov_wide$ELAC_L33=ELAC_cov_wide$Elac_L33_batch1+ELAC_cov_wide$Elac_L33_batch2
ELAC_cov_wide$ELAC_L15=ELAC_cov_wide$Elac_L15_batch1+ELAC_cov_wide$Elac_L15_batch2
ELAC_cov_wide$ELAC_L25=ELAC_cov_wide$Elac_L25_batch1+ELAC_cov_wide$Elac_L25_batch2
ELAC_cov_wide$ELAC_L35=ELAC_cov_wide$Elac_L35_batch1+ELAC_cov_wide$Elac_L35_batch2
colnames(ELAC_cov_wide)
for (i in 1:nrow(ELAC_cov_wide)){
  ELAC_cov_wide$ELAC3[i]=sum(as.numeric(ELAC_cov_wide[i,c(3,4,7,8,11,12)]),na.rm=TRUE)
  ELAC_cov_wide$ELAC5[i]=sum(as.numeric(ELAC_cov_wide[i,c(5,6,9,10,13,14)]),na.rm=TRUE)
}

# sum(as.numeric(ELAC_cov_wide[1,c(3,4,7,8,11,12)]),na.rm=TRUE)


WT_cov_wide <- spread(WT_cov, Type, score)
WT_cov_wide$WT_L13=WT_cov_wide$WT_L13_batch1+WT_cov_wide$WT_L13_batch2
WT_cov_wide$WT_L23=WT_cov_wide$WT_L23_batch1+WT_cov_wide$WT_L23_batch2
WT_cov_wide$WT_L33=WT_cov_wide$WT_L33_batch1+WT_cov_wide$WT_L33_batch2
WT_cov_wide$WT_L15=WT_cov_wide$WT_L15_batch1+WT_cov_wide$WT_L15_batch2
WT_cov_wide$WT_L25=WT_cov_wide$WT_L25_batch1+WT_cov_wide$WT_L25_batch2
WT_cov_wide$WT_L35=WT_cov_wide$WT_L35_batch1+WT_cov_wide$WT_L35_batch2
for (i in 1:nrow(WT_cov_wide)){
  WT_cov_wide$WT3[i]=sum(as.numeric(WT_cov_wide[i,c(3,4,7,8,11,12)]),na.rm=TRUE)
  WT_cov_wide$WT5[i]=sum(as.numeric(WT_cov_wide[i,c(5,6,9,10,13,14)]),na.rm=TRUE)
}


GFP_cov_wide <- spread(GFP_cov, Type, score)
GFP_cov_wide$GFP_L13=GFP_cov_wide$GFP_L13_batch1+GFP_cov_wide$GFP_L13_batch2
GFP_cov_wide$GFP_L23=GFP_cov_wide$GFP_L23_batch1+GFP_cov_wide$GFP_L23_batch2
GFP_cov_wide$GFP_L33=GFP_cov_wide$GFP_L33_batch1+GFP_cov_wide$GFP_L33_batch2
GFP_cov_wide$GFP_L15=GFP_cov_wide$GFP_L15_batch1+GFP_cov_wide$GFP_L15_batch2
GFP_cov_wide$GFP_L25=GFP_cov_wide$GFP_L25_batch1+GFP_cov_wide$GFP_L25_batch2
GFP_cov_wide$GFP_L35=GFP_cov_wide$GFP_L35_batch1+GFP_cov_wide$GFP_L35_batch2
for (i in 1:nrow(GFP_cov_wide)){
  GFP_cov_wide$GFP3[i]=sum(as.numeric(GFP_cov_wide[i,c(3,4,7,8,11,12)]),na.rm=TRUE)
  GFP_cov_wide$GFP5[i]=sum(as.numeric(GFP_cov_wide[i,c(5,6,9,10,13,14)]),na.rm=TRUE)
}

ELAC_long_no_repl <- gather(ELAC_cov_wide[,-c(3:20)], group, score, ELAC3:ELAC5, factor_key=TRUE)
ELAC_long_with_repl <- gather(ELAC_cov_wide[,-c(3:14,21:22)], group, score, ELAC_L13:ELAC_L35, factor_key=TRUE)
GFP_long_no_repl <- gather(GFP_cov_wide[,-c(3:20)], group, score, GFP3:GFP5, factor_key=TRUE)
GFP_long_with_repl <- gather(GFP_cov_wide[,-c(3:14,21:22)], group, score, GFP_L13:GFP_L35, factor_key=TRUE)
WT_long_no_repl <- gather(WT_cov_wide[,-c(3:20)], group, score, WT3:WT5, factor_key=TRUE)
WT_long_with_repl <- gather(WT_cov_wide[,-c(3:14,21:22)], group, score, WT_L13:WT_L35, factor_key=TRUE)

median(0,0,100)



new_all_cov_df_no_repl=rbind(ELAC_long_no_repl,GFP_long_no_repl,WT_long_no_repl)
new_all_cov_df_with_repl=rbind(ELAC_long_with_repl,GFP_long_with_repl,WT_long_with_repl)

new_all_cov_df_no_repl$group=as.character(new_all_cov_df_no_repl$group)
new_all_cov_df_no_repl$dpa=paste("dpa",substr(new_all_cov_df_no_repl$group,nchar(new_all_cov_df_no_repl$group),nchar(new_all_cov_df_no_repl$group)))
new_all_cov_df_no_repl$gene=substr(new_all_cov_df_no_repl$group,1,nchar(new_all_cov_df_no_repl$group)-1)



new_all_cov_df_no_repl %>%
  ggplot( aes(x=start, y=score, group=gene, color=gene)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE) +
  ggtitle("Normalized coverage") +
  theme_ipsum() +
  ylab("Normalized coverage")+
  facet_grid(~dpa)


all_high_cov_df=subset(new_all_cov_df_no_repl,new_all_cov_df_no_repl$score>2000000)
min(all_high_cov_df$start) #11760
max(all_high_cov_df$start) #13379
#peak is within 11760-13379
#12s rRNA location: 11674-12382
#16s rRNA location: 12590-13504

# new_all_cov_df_no_repl_part1=na.omit(subset(new_all_cov_df_no_repl,new_all_cov_df_no_repl$start<11674))
# new_all_cov_df_no_repl_part1=na.omit(subset(new_all_cov_df_no_repl,new_all_cov_df_no_repl$start<11674))
# new_all_cov_df_no_repl_12s=na.omit(subset(new_all_cov_df_no_repl,(new_all_cov_df_no_repl$start>=11674 & new_all_cov_df_no_repl$start<12383)))
# new_all_cov_df_no_repl_part2=na.omit(subset(new_all_cov_df_no_repl,(new_all_cov_df_no_repl$start>=12383 & new_all_cov_df_no_repl$start<12590)))
# new_all_cov_df_no_repl_16s=na.omit(subset(new_all_cov_df_no_repl,(new_all_cov_df_no_repl$start>=12590 & new_all_cov_df_no_repl$start<13505)))
# new_all_cov_df_no_repl_part3=na.omit(subset(new_all_cov_df_no_repl,new_all_cov_df_no_repl$start>13504))
# 
# #new_all_cov_df_no_repl_part1[which(new_all_cov_df_no_repl_part1$score==max(new_all_cov_df_no_repl_part1$score)),]
# max(new_all_cov_df_no_repl_part1$score)
# library(ggpubr)
# fragment1_no_repl_plot=new_all_cov_df_no_repl_part1 %>%
#   ggplot( aes(x=start, y=score, group=gene, color=gene)) +
#   geom_line() +
#   scale_color_viridis(discrete = TRUE) +
#   theme_ipsum() +
#   ylab("Normalized coverage")+
#   facet_grid(~dpa)
# fragment2_no_repl_plot=new_all_cov_df_no_repl_part2 %>%
#   ggplot( aes(x=start, y=score, group=gene, color=gene)) +
#   geom_line() +
#   scale_color_viridis(discrete = TRUE) +
#   theme_ipsum() +
#   ylab("Normalized coverage")+
#   facet_grid(~dpa)
# fragment3_no_repl_plot=new_all_cov_df_no_repl_part3 %>%
#   ggplot( aes(x=start, y=score, group=gene, color=gene)) +
#   geom_line() +
#   scale_color_viridis(discrete = TRUE) +
#   theme_ipsum() +
#   ylab("Normalized coverage")+
#   facet_grid(~dpa)
# 
# #?ggarrange
# ggarrange(fragment1_no_repl_plot,fragment2_no_repl_plot,fragment3_no_repl_plot,nrow=3,common.legend=TRUE)+ggtitle("Normalized coverage")
# 
# 
# all_cov_plot=new_all_cov_df_no_repl %>%
#   ggplot( aes(x=start, y=score, group=gene, color=gene)) +
#   geom_line() +
#   scale_color_viridis(discrete = TRUE) +
#   ggtitle("Normalized coverage") +
#   theme_ipsum() +
#   ylab("Normalized coverage")+
#   facet_grid(~dpa)

#mtDNA_anno_gtf=rtracklayer::import('E:/genome/mtDNA_anno_corrected.gtf')
#mtDNA_anno_gtf

library(rtracklayer)
library(GenomicFeatures)
#BiocManager::install("Gviz")
library(Gviz)
#BiocManager::install("txdbmaker")
library(txdbmaker)
# using the rtracklayer import function I can read the GFF3 from GenCode and create conversion table between gene id and gene symbol
#gff3 <- import.gff3("gencode.v24.annotation.gff3.gz")
gene2symbol <- mcols(mtDNA_anno_gtf)[,c("gene_id","gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id

# the GFF3 from GenCode can be parsed into TxDb object and after creating the GeneRegionTrack I could set the correct gene symbols using the table from previous step
txdb <- makeTxDbFromGFF('D:/genome/genome/mtDNA_anno_corrected.gtf')
geneTrack <- GeneRegionTrack(txdb,options(ucscChromosomeNames=FALSE))
?GeneRegionTrack
ranges(geneTrack)$symbol <- gene2symbol[ranges(geneTrack)$gene, "gene_name"]
ax <- GenomeAxisTrack()
# plotting with gene symbols
gen_REG_PLOT=plotTracks(list(ax, geneTrack),  showId=TRUE,geneSymbol=TRUE,
                        stacking="squish",collapseTranscripts=T)
?plotTracks
#plotTracks(txdb)


res <- plotTracks(list(ax, geneTrack), title.width = 0.6,  showId=TRUE)

#ggarrange(all_cov_plot,gen_REG_PLOT)

library(rtracklayer) 
#allChromosomeCoverage <- import.bw("data/small_Sorted_SRR568129.bw",as="GRanges") 
#ELAC_L131=import.bw("E:/Illumina/PARN_ELAC2_silencing/mRNA/mapped_to_chr_and_mtDNA/mtDNA_coverage/Elac_L13_batch1_coverage.txt",as="GRanges")
#----------------------------combine coverage and annotation--------------------
#part1: start=1, end=11673
#part2: start=12383, end=12589
#part3: start=13505, end=17175
#12s: start=11674, end=12382
#16s: start=12590, end=13504
# new_all_cov_df_no_repl
# all_no_re_plot_dpa3=new_all_cov_df_no_repl[new_all_cov_df_no_repl$dpa=="dpa 3",] %>%
#   ggplot( aes(x=start, y=score, group=gene, color=gene)) +
#   geom_line() +
#   scale_color_viridis(discrete = TRUE) +
#   theme_ipsum() +
#   ylab("Normalized coverage")+
#   ggtitle("Normalized coverage for dpa3")
# 
# 
# fragment1_no_repl_plot_dpa3=new_all_cov_df_no_repl_part1[new_all_cov_df_no_repl_part1$dpa=="dpa 3",] %>%
#   ggplot( aes(x=start, y=score, group=gene, color=gene)) +
#   geom_line() +
#   scale_color_viridis(discrete = TRUE) +
#   theme_ipsum() +
#   ylab("Normalized coverage")
# fragment2_no_repl_plot_dpa3=new_all_cov_df_no_repl_part2[new_all_cov_df_no_repl_part2$dpa=="dpa 3",] %>%
#   ggplot( aes(x=start, y=score, group=gene, color=gene)) +
#   geom_line() +
#   scale_color_viridis(discrete = TRUE) +
#   theme_ipsum() +
#   ylab("Normalized coverage")
# fragment3_no_repl_plot_dpa3=new_all_cov_df_no_repl_part3[new_all_cov_df_no_repl_part3$dpa=="dpa 3",] %>%
#   ggplot( aes(x=start, y=score, group=gene, color=gene)) +
#   geom_line() +
#   scale_color_viridis(discrete = TRUE) +
#   theme_ipsum() +
#   ylab("Normalized coverage")

library(cowplot)
library("grid")
library("ggplotify")
library(colorspace)
# 
# 
# res_part1 <- grid::grid.grabExpr(plotTracks(list(ax, geneTrack), title.width = 0.6,geneSymbol=TRUE, 
#                                             showId=TRUE, from=1, to=11673,collapseTranscripts=T))
# res_part1_gg=as.ggplot()
# res_part2 <- grid::grid.grabExpr(plotTracks(list(ax, geneTrack), title.width = 0.6, geneSymbol=TRUE,   
#                                             showId=TRUE,from=12383, to=12589,collapseTranscripts=T))
# res_part3 <- grid::grid.grabExpr(plotTracks(list(ax, geneTrack), title.width = 0.6,  geneSymbol=TRUE,  
#                                             showId=TRUE,from=13505, to=17175,collapseTranscripts=T))
# res_all <- grid::grid.grabExpr(plotTracks(list(ax, geneTrack), title.width = 0.05,geneSymbol=TRUE, 
#                                           showId=TRUE, collapseTranscripts=T))
# ?plotTracks

# 
# p1 = grid::grid.grabExpr(plotTracks(list(itrack, gtrack, grtrack), add = TRUE))
# p2 = ggplot2::qplot(1:10, 1:10)

# gridExtra::grid.arrange(fragment1_no_repl_plot_dpa3,res_part1,  ncol=1)
# gridExtra::grid.arrange(all_no_re_plot_dpa3,res_all,  ncol=1)
# plot_grid(all_no_re_plot_dpa3,res_all, ncol = 1, align = "h")
# ?grid.arrange
# 
# 
# grid.newpage()
# grid.draw(as.grob(all_no_re_plot_dpa3))
# vp = viewport(x=.35, y=.75, width=.35, height=.3)
# pushViewport(vp)
# grid.draw(res_all)
# upViewport()
# 
# 
# 
# 
# 
# grid.newpage()
# # 2x2 layout
# pushViewport(viewport(layout=grid.layout(2, 2)))
# # 1,1 first plot
# pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
# res_all
# #Gviz::plotTracks(list(ideogram, genomeaxis, siteA), add=TRUE)
# popViewport()
# # 2,2 second plot
# pushViewport(viewport(layout.pos.col=2,layout.pos.row=2))
# #Gviz::plotTracks(list(ideogram, genomeaxis, siteB), add=TRUE)
# popViewport()
# popViewport()

################################################################################
# #bigWig
# plotTracks(list(ax, geneTrack), title.width = 0.05,geneSymbol=TRUE, 
#            showId=TRUE, collapseTranscripts=T)
# allChromosomeCoverage <- import.bw("E:/Illumina/PARN_ELAC2_silencing/mRNA/mapped_to_chr_and_mtDNA/mtDNA_coverage_bigWig/Elac_L13_batch1_coverage.bw",as="GRanges")
# accDT <- DataTrack(allChromosomeCoverage,type="l") 
# plotTracks(list(ax,accDT,geneTrack),showId=TRUE) 
# ?DataTrack


options(max.print=17000)
getOption("max.print")
options(max.print = .Machine$integer.max)
#convert df into GRanges:
ELAC_long_no_repl_dpa3=subset(ELAC_long_no_repl,ELAC_long_no_repl$group=="ELAC3")
WT_long_no_repl_dpa3=subset(WT_long_no_repl,WT_long_no_repl$group=="WT3")
GFP_long_no_repl_dpa3=subset(GFP_long_no_repl,GFP_long_no_repl$group=="GFP3")
colnames(ELAC_long_no_repl_dpa3)=c("seqnames","start","group","ELAC3_cov")
colnames(WT_long_no_repl_dpa3)=c("seqnames","start","group","WT3_cov")
colnames(GFP_long_no_repl_dpa3)=c("seqnames","start","group","GFP3_cov")
dpa3_merged_no_repl=merge(ELAC_long_no_repl_dpa3[,-3],WT_long_no_repl_dpa3[,-3],by=c("start","seqnames"))
dpa3_merged_no_repl=merge(dpa3_merged_no_repl,GFP_long_no_repl_dpa3[,-3],by=c("start","seqnames"))
dpa3_merged_no_repl$end=dpa3_merged_no_repl$start+1
dpa3_merged_no_repl=dpa3_merged_no_repl[order(dpa3_merged_no_repl$start),]
dpa3_merged_no_repl_GR=makeGRangesFromDataFrame(dpa3_merged_no_repl,keep.extra.columns=TRUE)


ELAC_long_no_repl_dpa5=subset(ELAC_long_no_repl,ELAC_long_no_repl$group=="ELAC5")
WT_long_no_repl_dpa5=subset(WT_long_no_repl,WT_long_no_repl$group=="WT5")
GFP_long_no_repl_dpa5=subset(GFP_long_no_repl,GFP_long_no_repl$group=="GFP5")
colnames(ELAC_long_no_repl_dpa5)=c("seqnames","start","group","ELAC5_cov")
colnames(WT_long_no_repl_dpa5)=c("seqnames","start","group","WT5_cov")
colnames(GFP_long_no_repl_dpa5)=c("seqnames","start","group","GFP5_cov")
dpa5_merged_no_repl=merge(ELAC_long_no_repl_dpa5[,-3],WT_long_no_repl_dpa5[,-3],by=c("start","seqnames"))
dpa5_merged_no_repl=merge(dpa5_merged_no_repl,GFP_long_no_repl_dpa5[,-3],by=c("start","seqnames"))
dpa5_merged_no_repl$end=dpa5_merged_no_repl$start+1
dpa5_merged_no_repl=dpa5_merged_no_repl[order(dpa5_merged_no_repl$start),]
dpa5_merged_no_repl_GR=makeGRangesFromDataFrame(dpa5_merged_no_repl,keep.extra.columns=TRUE)

ax <- GenomeAxisTrack()


subset(ELAC_long_no_repl_dpa3,ELAC_long_no_repl_dpa3$start>7280 & ELAC_long_no_repl_dpa3$start<7300)
?DataTrack
accDT3 <- DataTrack(dpa3_merged_no_repl_GR,type="a",groups=c("ELAC_dpa3","WT_dpa3","GFP_dpa3")) 
accDT5 <- DataTrack(dpa5_merged_no_repl_GR,type="a",groups=c("ELAC_dpa5","WT_dpa5","GFP_dpa5")) 
#group(accDT)=c("ELAC_dpa3","WT_dpa3","GFP_dpa3")
#?DataTrack

# 
# 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE) 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=1, to=11673) 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=12383, to=12589) 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=13505, to=17175) 



options(showTailLines=Inf)
mtDNA_anno_gtf

#select groups
#gr1: from=2250, to=2450
#gr2: from=4440, to=4700
#gr3: from=5450, to=5650
#gr4: from=6250, to=6350
#gr5: from=

# 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=2250, to=2450) 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=4440, to=4700) 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=6250, to=6350) 


?plotTracks
?HighlightTrack
# ht1 <- HighlightTrack(trackList = accDT3, start = c(2200, 4440),end=c(2450,4700))
# ht2 <- HighlightTrack(trackList = accDT5, start = c(2200, 4440),end=c(2450,4700))
# ht3 <- HighlightTrack(trackList = geneTrack, start = c(2200, 4440),end=c(2450,4700))

# plotTracks(list(ax,ht1,ht2,ht3),showId=TRUE,just.group = "above",
#            CDS = "darkred", ncRNA = "darkgreen")
names(accDT3)="Normalized coverage"
names(accDT5)="Normalized coverage"
#reg1
?plotTracks
reg1=plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=2350, to=2500,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "above",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 1",cex.main = 3) 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=2370, to=2420,fontcolor.feature = "darkblue",
#            name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "above",
#            CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 1",cex.main = 3) 
#reg2
plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=4400, to=4700,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 2",cex.main = 3) 
reg2=plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=4600, to=4700,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 2",cex.main = 3) 
#reg3
#Figure_5_2
plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=5600, to=5700,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 3",cex.main = 3) 
#reg4
#Figure_5_3
plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=6250, to=6450,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 4",cex.main = 3) 
plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=6300, to=6450,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 4",cex.main = 3) 
#reg5
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=7250, to=7400,fontcolor.feature = "darkblue",
#            name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
#            CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 5",cex.main = 3) 
plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=7280, to=7330,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 5",cex.main = 3) 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=7280, to=7320,fontcolor.feature = "darkblue",
#            name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
#            CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 5",cex.main = 3) 
#reg6
plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=8300, to=8400,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 6",cex.main = 3)
#reg7
plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=8750, to=8850,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 7",cex.main = 3) 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=8760, to=8790,fontcolor.feature = "darkblue",
#            name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
#            CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 7",cex.main = 3) 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=11650, to=11700,fontcolor.feature = "darkblue",
#            name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
#            CDS = "darkred", ncRNA = "darkgreen") 
#reg8
#Figure_5_3
plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=11470, to=11685,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 8",cex.main = 3) 
#reg8_1
plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=11470, to=12200,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 8_1",cex.main = 3) 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=12380, to=13000,fontcolor.feature = "darkblue",
#            name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "above",
#            CDS = "darkred", ncRNA = "darkgreen")
#reg9
#Figure_5_9
plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=13505, to=13800,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 9",cex.main = 3) 
plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=13700, to=13800,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 9",cex.main = 3) 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=13505, to=15000,fontcolor.feature = "darkblue",
#            name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
#            CDS = "darkred", ncRNA = "darkgreen") 
#reg10
plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=14450, to=14600,fontcolor.feature = "darkblue",
           name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
           CDS = "darkred", ncRNA = "darkgreen",showFeatureId=TRUE,main="Region 10",cex.main = 3) 
# plotTracks(list(ax,accDT3,accDT5,geneTrack),showId=TRUE, from=15550, to=16000,fontcolor.feature = "darkblue",
#            name = "Normalized coverage",featureAnnotation = "gene_name",just.group = "left",
#            CDS = "darkred", ncRNA = "darkgreen") 





?HighlightTrack

ht <- HighlightTrack(trackList = c(accDT3,accDT5,geneTrack), start = c(2350, 4400,5600,6250,7250,8300,8750,11470,13505,14450),
                     end=c(2500,4700,5700,6450,7400,8400,8850,12200,13800,14600))
#initialize(ht)
#length(ht)
#?plotTracks
plotTracks(list(ax,ht),showId=TRUE,just.group = "above",
           CDS = "darkred", ncRNA = "darkgreen")
?just.group

feature(geneTrack)

#EXAMPLE
################################################################################
## Create some tracks to plot
st <- c(2000000, 2070000, 2100000, 2160000)
ed <- c(2050000, 2130000, 2150000, 2170000)
str <- c("-", "+", "-", "-")
gr <- c("Group1", "Group2", "Group1", "Group3")
annTrack <- AnnotationTrack(
  start = st, end = ed, strand = str, chromosome = 7,
  genome = "hg19", feature = "test", group = gr,
  id = paste("annTrack item", 1:4),
  name = "annotation track foo",
  stacking = "squish"
)

ax <- GenomeAxisTrack()

# dt <- DataTrack(
#   start = seq(min(st), max(ed), len = 10), width = 18000,
#   data = matrix(runif(40), nrow = 4), genome = "hg19", chromosome = 7,
#   type = "histogram", name = "data track bar"
# )
# 
# 
# 
# ## Now plot the tracks
# res <- plotTracks(list(ax, annTrack, dt),  showId=TRUE)
# 
# ## Plot only a subrange
# res <- plotTracks(list(ax, annTrack, dt), from = 2080000, to = 2156000)
# 
# ## Extend plotting ranges
# res <- plotTracks(list(ax, annTrack, dt), extend.left = 200000, extend.right = 200000)
# 
# ## Add a header
# res <- plotTracks(list(ax, annTrack, dt),
#                   main = "A GenomGraphs plot",
#                   col.main = "darkgray"
# )
# 
# ## Change vertical size and title width
# res <- plotTracks(list(ax, annTrack, dt), sizes = c(1, 1, 5))

names(annTrack) <- "foo"
res <- plotTracks(list(ax, annTrack), title.width = 0.6,  showId=TRUE)

## Adding and lattice like plots
library(grid)
grid.newpage()
pushViewport(viewport(height = 0.5, y = 1, just = "top"))
grid.rect()
plotTracks(annTrack, add = TRUE)
popViewport(1)
pushViewport(viewport(height = 0.5, y = 0, just = "bottom"))
grid.rect()
plotTracks(dt, add = TRUE)
popViewport(1)
## Not run: 
library(lattice)
myPanel <- function(x, ...) {
  plotTracks(annTrack,
             panel.only = TRUE,
             from = min(x), to = max(x), shape = "box"
  )
}
a <- seq(1900000, 2250000, len = 40)
xyplot(b ~ a | c, data.frame(a = a, b = 1, c = cut(a, 4)),
       panel = myPanel,
       scales = list(x = "free")
)
################################################################################




library(ggtranscript)

sod1_annotation %>% head()
# extract exons
sod1_exons <- sod1_annotation %>% dplyr::filter(type == "exon")
sod1_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = transcript_biotype)
  ) +
  geom_intron(
    data = to_intron(sod1_exons, "transcript_name"),
    aes(strand = strand)
  )

sod1_rescaled <- shorten_gaps(
  sod1_exons, 
  to_intron(sod1_exons, "transcript_name"), 
  group_var = "transcript_name"
)

sod1_rescaled %>%
  dplyr::filter(type == "exon") %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = transcript_biotype)
  ) +
  geom_intron(
    data = sod1_rescaled %>% dplyr::filter(type == "intron"), 
    arrow.min.intron.length = 200
  )

sod1_annotation %>% head()
# extract exons
sod1_exons <- sod1_annotation %>% dplyr::filter(type == "exon")
sod1_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = transcript_biotype)
  ) +
  geom_intron(
    data = to_intron(sod1_exons, "transcript_name"),
    aes(strand = strand)
  )

mtDNA_anno_gtf_tibbble=as_tibble(mtDNA_anno_gtf)
#mtDNA_anno_gtf_tibbble$type
mtDNA_anno_gtf_tibbble %>% head()
# extract exons
mtDNA_anno_gtf_tibbble <- mtDNA_anno_gtf_tibbble %>% dplyr::filter(type == "transcript")
mtDNA_anno_gtf_tibbble %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = gene_name
  )) +
  geom_intron(
    data = to_intron(mtDNA_anno_gtf_tibbble, "gene_name"),
    aes(strand = strand)
  )













library("rtracklayer")
library("ggcoverage")

meta.file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
sample.meta = read.csv(meta.file)
# track folder
track.folder = system.file("extdata", "RNA-seq", package = "ggcoverage")
# load bigwig file
track.df = LoadTrackFile(track.folder = track.folder, format = "bw",
                         meta.info = sample.meta)
# check data
head(track.df)
tail(track.df)
gtf.file = system.file("extdata", "used_hg19.gtf", package = "ggcoverage")
gtf.gr = rtracklayer::import.gff(con = gtf.file, format = 'gtf')
?ggcoverage
(basic.coverage = ggcoverage(data = track.df, color = "auto", 
                             range.position = "out")+geom_gene(gtf.gr=gtf.gr))




head(new_all_cov_df_no_repl_part1)
head(track.df)
new_all_cov_df_no_repl_part1_ggcov=new_all_cov_df_no_repl_part1
new_all_cov_df_no_repl_part1_ggcov$end=new_all_cov_df_no_repl_part1_ggcov$start+1
colnames(new_all_cov_df_no_repl_part1_ggcov)=c("seqnames","start","Type" ,"score", "dpa","Group", "end" )
ggcoverage(data = new_all_cov_df_no_repl_part1_ggcov, color = "auto", 
           range.position = "out")+geom_gene(gtf.gr=mtDNA_anno_gtf)
?geom_gene
mtDNA_anno_gtf

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2")
install.packages("TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggbio)
data(genesymbol, package = "biovizBase")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

autoplot(mtDNA_anno_gtf, aes(type = gene_name))
autoplot(mtDNA_anno_gtf)+ theme_genome()
p.txdb
?autoplot


p17 <- autoplot(mtDNA_anno_gtf)
plotRangesLinkedToData(mtDNA_anno_gtf, annotation = list(p17))

?plotRangesLinkedToData

exon.new=mtDNA_anno_gtf
values(exon.new)$sample1 <- rnorm(length(exon.new), 10, 3)
values(exon.new)$sample2 <- rnorm(length(exon.new), 10, 10)
values(exon.new)$score <- rnorm(length(exon.new))
p17 <- autoplot(exon.new)
plotRangesLinkedToData(mtDNA_anno_gtf,stat.y = c("sample1", "sample2"),annotation = list(p17))








# Plot
# ELAC_WT %>%
#   ggplot( aes(x=start, y=score, group=Type, color=Type)) +
#   geom_line() +
#   scale_color_viridis(discrete = TRUE) +
#   ggtitle("Normalized coverage") +
#   theme_ipsum() +
#   ylab("Normalized coverage")
all_cov_df$set=substr(all_cov_df$Type,1,nchar(all_cov_df$Type)-7)
all_cov_df$dpa=paste("dpa",substr(all_cov_df$Type,nchar(all_cov_df$set),nchar(all_cov_df$set)))
all_cov_df$gene=substr(all_cov_df$set,1,nchar(all_cov_df$set)-4)
all_cov_df %>%
  ggplot( aes(x=start, y=score, group=gene, color=gene)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE) +
  ggtitle("Normalized coverage") +
  theme_ipsum() +
  ylab("Normalized coverage")+
  facet_grid(~dpa)

all_low_cov_df=subset(all_cov_df,all_cov_df$score<2000000)
all_high_cov_df=subset(all_cov_df,all_cov_df$score>2000000)


all_high_cov_df %>%
  ggplot( aes(x=start, y=score, group=gene, color=gene)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE) +
  ggtitle("Normalized coverage") +
  theme_ipsum() +
  ylab("Normalized coverage")+
  facet_grid(~dpa)




tail(ELAC_cov_wide$start)
histogram(ELAC_cov_wide$start)
install.packages("plotgardener")
library(plotgardener)




library(ggbio)

# load gene annotations in GTF format
genes <- mtDNA_anno_gtf

# create a GRanges object with gene annotations
library(GenomicRanges)
gr <- makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)

# load coverage data for each sample in BED format
cov1 <- read.table("sample1.bed")
cov2 <- read.table("sample2.bed")
# ...
# create Coverage objects for each sample
?coverage
covobj1 <- coverage(ELAC_L131, gr)
covobj2 <- coverage(cov2, gr)
# ...

# merge gene annotations with coverage data for each sample
covgr1 <- merge(covobj1, gr)
covgr2 <- merge(covobj2, gr)
# ...

# compute normalized coverage for each gene
covnorm1 <- as.data.frame(normalize(covgr1))
covnorm2 <- as.data.frame(normalize(covgr2))
# ...
# create a ggplot object for each sample
gg1 <- ggplot(covnorm1, aes(x = start, y = name)) +
  geom_rect(aes(xmin = start, xmax = end, ymin = name - 0.4, ymax = name + 0.4,
                fill = score), color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  xlab("Position") + ylab("Gene") + ggtitle("Sample 1")

gg2 <- ggplot(covnorm2, aes(x = start, y = name)) +
  geom_rect(aes(xmin = start, xmax = end, ymin = name - 0.4, ymax = name + 0.4,
                fill = score), color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  xlab("Position") + ylab("Gene") + ggtitle("Sample 2")
# ...

# arrange the ggplot objects in a grid
library(cowplot)
plot_grid(gg1, gg2, ..., ncol = 1)



#install.packages("Gviz")
library(Gviz)

# Load gene annotation
gtf <- mtDNA_anno_gtf

# Filter for exons
exons <- gtf[gtf$V3 == "exon",]

# Create a genome axis track
genome_track <- GenomeAxisTrack()
strand(gtf)
# Create an annotation track
?AnnotationTrack
gene_track <- AnnotationTrack(range = gtf[,c(1,4,5)],
                              feature = gtf$gene_name,
                              name = "Genes",
                              featureType = "transcript",
                              strand = strand(gtf),
                              fill = "#FFC0CB",
                              stack = TRUE,
                              options(ucscChromosomeNames=FALSE))

# Create a coverage track
cov_track <- DataTrack(range = exons[,c(1,4,5)],
                       data = list(file1, file2, file3, ...),
                       name = c("Sample1", "Sample2", "Sample3", ...),
                       type = "l",
                       fill = c("black", "red", "blue", ...),
                       ylim = c(0, max_coverage),
                       from = min_position,
                       to = max_position)

# Plot tracks
plotTracks(list(genome_track, gene_track, cov_track))






## Coverage of a GRanges object:
gr <- GRanges(
  seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges=IRanges(1:10, end=10),
  strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  seqlengths=c(chr1=11, chr2=12, chr3=13))
cvg <- coverage(gr)
pcvg <- coverage(gr[strand(gr) == "+"])
mcvg <- coverage(gr[strand(gr) == "-"])
scvg <- coverage(gr[strand(gr) == "*"])
stopifnot(identical(pcvg + mcvg + scvg, cvg))

## Coverage of a GPos object:
pos_runs <- GRanges(c("chr1", "chr1", "chr2"),
                    IRanges(c(1, 5, 9), c(10, 8, 15)))
gpos <- GPos(pos_runs)
coverage(gpos)

## Coverage of a GRangesList object:
gr1 <- GRanges(seqnames="chr2",
               ranges=IRanges(3, 6),
               strand = "+")
gr2 <- GRanges(seqnames=c("chr1", "chr1"),
               ranges=IRanges(c(7,13), width=3),
               strand=c("+", "-"))
gr3 <- GRanges(seqnames=c("chr1", "chr2"),
               ranges=IRanges(c(1, 4), c(3, 9)),
               strand=c("-", "-"))
grl <- GRangesList(gr1=gr1, gr2=gr2, gr3=gr3)



