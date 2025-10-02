
# planarian-elac2-scripts

This repository contains all scripts used in the analysis for:

**Knockdown of Smed ELAC2 in *Schmidtea mediterranea* results in delayed regeneration and a reduced accumulation of the small non-coding RNA, 5′ tiRNA-Gly-GCC**  

Scripts are organized by figure (main and supplementary).

---

## Main Figures

### Figure 1 – Annotation of Smed ELAC2
- `fig1_annotation.R` – Assembly of Smed ELAC2 transcript, comparison with PlanMine models, IGV plots.  
- `sequence_analysis.R` – Protein domain and mitochondrial targeting predictions.

### Figure 2 – Eye regeneration phenotypes
- `fig2_eye_regeneration.R` – qPCR knockdown validation, quantification of eye defects, plotting regeneration delay.

### Figure 3 – Nerve cord regeneration
- `fig3_nerve_regeneration.R` – FISH quantification of PC2 marker, image processing and fluorescence intensity analysis.

### Figure 4 – Cell population stability
- `fig4_flow_cytometry.R` – Flow cytometry gating, quantification of X1/X2/Xins cell populations.

### Figure 5 – Mitochondrial transcript processing
- `fig5_mt_transcripts.R` – Coverage of mitochondrial junction regions, qPCR validation of unprocessed transcripts.  
- `mtDNA_chimeras.R` – PCR validation of mtDNA chimeras (junction assays).

### Figure 6 – Gene expression profiling
- `fig6_gene_expression.R` – DESeq2 analysis of long RNA-seq, GO term enrichment, cell-type assignment.

### Figure 7 – Small RNA landscape
- `fig7_sncRNA_landscape.R` – Distribution of small RNA classes, differential accumulation analysis.  
- `sncRNA_analysis.R` – Core small RNA processing pipeline (adapter trimming, mapping, classification).

### Figure 8 – ER stress analysis
- `fig8_er_stress.R` – Heatmap generation of ER stress–related genes.  
- `ER_stress.R` – Alternative pipeline for expression of IRE1α pathway genes.

### Figure 9 – 5′ tiRNA-Gly-GCC functional validation
- `fig9_tirna_target.R` – Target validation plots for protein tyrosine phosphatase (PTP).  
- `soaking_experiment.R` – Analysis of rescue experiments with synthetic 5′ tiRNA-Gly-GCC mimic.  
- `chimeras_PCR.R` – PCR validation of predicted miRNA-like interactions.

---

## Supplementary Figures and Data

- `supplement_mtDNA_validation.R` – Extended analysis of unprocessed mt-tRNA junctions.  
- `supplement_chimera_analysis.R` – In silico validation of chimeric reads.  
- `supplement_GO_analysis.R` – Extended GO clustering for DEGs.  
- `supplement_primer_sets.R` – Primer design and validation checks for qRT-PCR and FISH probes.  

---

## General Pipelines

- `pipeline_longRNAseq.R` – Long RNA-seq pipeline (QC, STAR mapping, StringTie quantification).  
- `pipeline_smallRNAseq.R` – Small RNA-seq pipeline (QC, trimming, Bowtie2 mapping).  
- `DESeq2_analysis.R` – Differential gene and small RNA expression analysis.  
- `GO_analysis.R` – GO enrichment and clustering using topGO.  
- `miRanda_prediction.py` – miRNA-like target prediction of 5′ tiRNA-Gly-GCC.  
- `ufold_prediction.py` – Secondary structure filtering for candidate targets.  

---

## Data Availability

RNA-seq data: [GSE293999](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293999), [GSE294003](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294003), [GSE294000](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294000) (BioProject [PRJNA1243063](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1243063))  
Manuscript DOI: (to be added upon acceptance)

---

## Citation

If you use these scripts, please cite:  
> Dutta A., Zaremba A., Samelak-Czajka A., Marszalek-Zenczak M., Trybus M., Osuch M., Figlerowicz M., Jackowiak P. Knockdown of *Smed ELAC2* in *Schmidtea mediterranea* results in delayed regeneration and a reduced accumulation of the small non-coding RNA, 5′ tiRNA-Gly-GCC. Nucleic Acids Research (2025).

---

