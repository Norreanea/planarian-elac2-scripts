############################################################
# ELAC2 gene model visualization 
############################################################

# ---- Packages ----
# If needed, install Bioc packages once:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GenomicFeatures", "rtracklayer", "Gviz"))

library(GenomicFeatures)  
library(rtracklayer)      
library(Gviz)             

# ---- Point to the GTF  ----
gtf_url <- "https://raw.githubusercontent.com/Norreanea/planarian-elac2-scripts/main/data/ELAC2_structure_comparison.gtf "

# ---- Download to a temp file  ----
gtf_tmp <- tempfile(fileext = ".gtf")
download.file(gtf_url, destfile = gtf_tmp, mode = "wb")

# ---- Build TxDb and import annotation ----
txdb <- makeTxDbFromGFF(gtf_tmp)
gtf_gr <- rtracklayer::import(gtf_tmp)  

# ---- Build Gviz tracks ----
# Turn off UCSC-style chromosome name enforcement 
options(ucscChromosomeNames = FALSE)

geneTrack <- GeneRegionTrack(txdb,
                             name   = "ELAC2 locus",
                             showId = TRUE)

axisTrack <- GenomeAxisTrack()

# Tweak visual parameters on the gene track
displayPars(geneTrack)$cex.group        <- 2    # group label size
displayPars(geneTrack)$groupAnnotation  <- "group"

# ---- Plot ----
# pdf("elac2_gene_model.pdf", width = 8, height = 4)
plotTracks(
  list(axisTrack, geneTrack),
  showId   = TRUE,
  cex      = 2,    # overall text size
  cex.axis = 3,    # axis labels
  cex.main = 7,   
  cex.title= 3     # track title size
)
# dev.off()
