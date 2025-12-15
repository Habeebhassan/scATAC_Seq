# scATAC-seq Analysis Project - SCB 2025-26

**Course:** Single Cell Bioinformatics 2025-26 (Müller Lab @ Saarland University)  
**Project:** Project 2 - Single Cell ATAC-seq Analysis w/ ArchR 
**Genome Assembly:** hg38 (GRCh38)

## Project Overview

This repository contains the analysis pipeline for single-cell ATAC-seq (scATAC-seq) data using the **ArchR** framework in R. The goal of the project is to process raw fragment files, assess quality, remove batch effects, and identify cell-type-specific chromatin accessibility profiles.

The analysis is based on three biological replicates:
* `ATAC_555_2`
* `ATAC_557`
* `ATAC_HIP02_frozen`

## Pipeline & Methodology

The analysis follows a structured weekly workflow:

### [cite_start]Week 1: Preprocessing & Quality Control [cite: 26]
* [cite_start]**Data Import:** Raw fragment files were mapped to the **hg38** reference genome[cite: 30].
* **Quality Control (QC):**
  * [cite_start]**Fragment Size Distribution:** Validated nucleosomal periodicity (<100bp, ~200bp, ~400bp)[cite: 47].
  * [cite_start]**TSS Enrichment:** Assessed signal-to-noise ratio (High quality: Score > 8)[cite: 50].
* [cite_start]**Filtering:** Applied strict thresholds (TSS > 8, Fragments > 1000) to remove low-quality cells[cite: 60].
* [cite_start]**Doublet Removal:** Identified and removed heterotypic doublets using ArchR's simulation method[cite: 38].

### [cite_start]Week 2: Dimensionality Reduction & Clustering [cite: 65]
* [cite_start]**Peak Calling:** Generated a reproducible peak set using **MACS2** on pseudo-bulked replicates grouped by sample[cite: 67, 71].
* [cite_start]**Dimensionality Reduction:** Performed Iterative Latent Semantic Indexing (**LSI**) on the peak matrix, excluding the first component (correlation with depth)[cite: 79].
* [cite_start]**Batch Correction:** Corrected significant sample-specific batch effects using **Harmony**[cite: 88].
* [cite_start]**Clustering:** Applied **Louvain clustering** on the Harmony-corrected space, resulting in 5 distinct, mixed clusters[cite: 91].
* [cite_start]**Marker Discovery:** Identified cluster-specific accessible peaks to define regulatory signatures[cite: 74].

## Repository Structure
