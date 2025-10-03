
# planarian-elac2-scripts

This repository contains scripts used in the analysis for:

**Knockdown of Smed ELAC2 in *Schmidtea mediterranea* results in delayed regeneration and a reduced accumulation of the small non-coding RNA, 5′ tiRNA-Gly-GCC**  

Scripts are organized by figure (main and supplementary).

---

## Main Figures

### Figure 1 – Annotation of Smed ELAC2
- `scripts/main/fig1_annotation.R` – Comparison with PlanMine models, IGV-style visualization.  

### Figure 2 – Eye regeneration phenotypes
- `scripts/main/fig2_phenotype.R` – qPCR knockdown validation, plotting regeneration delay.

### Figure 5 – Mitochondrial transcript processing
- `scripts/main/fig5_mtDNA_coverage.R` – Coverage of mitochondrial junction regions.  

### Figure 6 – Gene expression profiling
- `scripts/main/fig6_DEG.R` – DESeq2 analysis of long RNA-seq, GO term enrichment, cell-type assignment.

### Figure 7 – Small RNA landscape
- `scripts/main/fig7_snRNA_analysis.R` – Distribution of small RNA classes, differential accumulation analysis.  

### Figure 8 – ER stress analysis
- `scripts/main/fig8_ER_stress.R` – Heatmap generation of ER stress–related genes.  

### Figure 9 – Rescue experiments with 5′ tiRNA-Gly-GCC
- `scripts/main/fig9_soaking.R` – Target validation plots for protein tyrosine phosphatase (PTP).
- [miRanda_prediction](https://github.com/Norreanea/miRNA-seed-matching-counter-for-miRanda-output) - miRNA-like target prediction.

---

## Supplementary Figures and Data

- `scripts/supplementary/figS5_ncRNA_length_and_tRF_scheme.R` – Generates length-distribution plots for ncRNA classes and a schematic of tRNA-derived fragment.  
  

---


## Data Availability

RNA-seq data: [GSE293999](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293999), [GSE294003](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294003), [GSE294000](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294000) (BioProject [PRJNA1243063](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1243063))  


---

## Citation

If you use these scripts, please cite:  
> Dutta A., Zaremba A., Samelak-Czajka A., Marszalek-Zenczak M., Trybus M., Osuch M., Figlerowicz M., Jackowiak P. Knockdown of *Smed ELAC2* in *Schmidtea mediterranea* results in delayed regeneration and a reduced accumulation of the small non-coding RNA, 5′ tiRNA-Gly-GCC. Manuscript DOI: (to be added upon acceptance)

---

