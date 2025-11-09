# Differential Gene Expression Analysis in Lewy Body Disease (LBD)

**Project:** PRJNA1023207
**Organism:** *Homo sapiens*
**Tissue:** Anterior Cingulate Cortex
**Sequencing Strategy:** Bulk RNA-Seq
**Analysis Toolchain:** R / DESeq2 / enrichplot

---

## Overview

Lewy Body Disease (LBD) and Alzheimer’s Disease (AD) are major neurodegenerative disorders that contribute to cognitive impairment.

* **LBD** is characterized by the accumulation of **α-synuclein** in Lewy bodies and neurites, often co-occurring with tau tangle pathology.
* **AD** is primarily associated with **amyloid plaques** and **neurofibrillary tangles**.

Both diseases can present overlapping pathologies, including vascular and TDP-43 changes. This project investigates **transcriptional changes in the anterior cingulate cortex**, a brain region involved in cognitive control and emotional regulation, to identify shared and distinct molecular features in LBD.

---

## Dataset Summary

This study uses a subset of RNA-Seq samples from **BioProject [PRJNA1023207](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1023207)**.
All samples were obtained from the **anterior cingulate cortex** and sequenced using the **Illumina platform**.

### Sample Overview

| Group   | Sex    | Sample Count | Example Run Accessions                                          |
| ------- | ------ | ------------ | --------------------------------------------------------------- |
| LBD     | Female | 5            | SRR26254070, SRR26254010, SRR26254144, SRR26254284, SRR26253954 |
| LBD     | Male   | 5            | SRR26254115, SRR26254278, SRR26253870, SRR26254409, SRR26253979 |
| Control | Female | 2            | SRR26254419, SRR26254262                                        |
| Control | Male   | 2            | SRR26254102, SRR26254264                                        |

**Total Samples:** 14
**Library Type:** RNA-Seq (paired-end)
**Organism:** *Homo sapiens*

---

## Analysis Workflow

### Data Preprocessing

* Raw read counts imported from `featureCounts` output (`counts_clean.tsv`).
* Metadata table (`samples_metadata.tsv`) defines:
  `sample`, `disease_status (LBD/Normal)`, `sex (Male/Female)`, and sequencing run.
* Filtered genes with low counts: `rowSums(counts >= 10) >= 2`.

### Differential Expression (DE) Analysis

* Conducted using **DESeq2** with the model:

  ```r
  design = ~ sex + disease_status
  ```
* Comparison: `LBD vs Normal`
* Result: 19,284 genes tested.
* Identified **~16 DEGs** at FDR < 0.1 (7 upregulated, 9 downregulated).

### Normalization & Visualization

* Variance Stabilizing Transformation (`vst`) for normalization.
* **Visualizations Generated:**

  * Volcano plot (log2FC vs -log10 padj)
  * PCA plot (sex and disease status separation)
  * Heatmap of top 30 DEGs (z-score normalized)

### Functional Enrichment Analysis

Performed using **clusterProfiler** and **org.Hs.eg.db**.

* **Gene Ontology (GO Biological Process):**

  * Enriched in processes such as *positive regulation of vasoconstriction*, *maintenance of cell polarity*, and *amino/aminoglycan catabolic process*.
* **KEGG Pathways:**

  * Significant enrichment in *Mucin type O-glycan biosynthesis*, *Oxytocin signaling*, *Calcium/cAMP signaling*, and *Neuroactive ligand-receptor interaction*.

---

## Key Visualizations

| Plot             | Description                                                                                |
| ---------------- | ------------------------------------------------------------------------------------------ |
| **Volcano Plot** | Highlights up/downregulated genes between LBD and controls.                                |
| **PCA Plot**     | Displays partial separation between LBD and control samples; visible sex-related variance. |
| **Heatmap**      | Top 30 DEGs clustered by expression z-scores, showing disease-driven expression patterns.  |
| **GO Dotplot**   | Biological processes enriched among significant DEGs.                                      |
| **KEGG Dotplot** | Pathway-level functional enrichment visualization.                                         |

---

## Findings Summary

* The differential expression signal between LBD and control samples is modest, suggesting subtle transcriptional alterations.
* Despite few DEGs, pathway analysis reveals **consistent enrichment in signaling and glycosylation processes**, both relevant to neurodegenerative disease mechanisms.
* Sex adjustment partially clarified clustering, indicating minor sex-associated transcriptomic differences.

---

## Tools & Packages

| Tool / Library         | Purpose                               |
| ---------------------- | ------------------------------------- |
| **DESeq2**             | Differential gene expression analysis |
| **clusterProfiler**    | GO and KEGG enrichment analysis       |
| **org.Hs.eg.db**       | Gene annotation mapping               |
| **enrichplot**         | Visualization of enrichment results   |
| **ggplot2 / pheatmap** | Data visualization                    |

---
### How to Run This Project

```bash
# Clone the repository
git clone https://github.com/<yourusername>/RNA_SEQ.git
cd RNA_SEQ

# Create and activate the conda environment
conda env create -f environment.yml
conda activate rnaseq_env

# Run the differential expression analysis in R
Rscript scripts/LBD_vs_Normal.R

# (Optional) Run the Nextflow pipeline for alignment and counting
nextflow run scripts/rna_seq.nf -profile conda
```


## Output Files

| File                                          | Description                                |
| --------------------------------------------- | ------------------------------------------ |
| `DESeq2_LBD_vs_Normal_results.csv`            | Raw DESeq2 results table                   |
| `GO_BP_enrich_results_relaxed.csv`            | GO enrichment results (Biological Process) |
| `KEGG_LBD_vs_Normal.png`                      | KEGG enrichment dotplot                    |
| `GO_LBD_vs_Normal.png`                        | GO enrichment dotplot                      |
| `Volcano_plot.png`                            | Volcano plot of DEGs                       |
| `pca_LBD_vs_Normal_sexAdjusted.png`           | PCA visualization                          |
| `heatmap_top30_LBD_vs_Normal_sexAdjusted.png` | Top 30 DEGs heatmap                        |

---
