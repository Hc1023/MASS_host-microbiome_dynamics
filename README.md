# MASS_host-microbiome_dynamics

This repository contains the analysis code and processed outputs for the study:

**“Longitudinal host-microbiome dynamics define immune trajectory divergence and mortality risk in severe pneumonia”**

## Overview

Severe pneumonia often progresses to sepsis and is associated with high mortality in ICU settings. This project integrates **whole-blood metatranscriptomics** and **clinical phenotyping** to characterize longitudinal host-microbe interactions.

Using a prospective cohort of **417 patients and 976 samples**, we:

* Identified **divergent host transcriptional trajectories**
* Characterized **circulating microbial dynamics**
* Linked **herpesvirus reactivation (HCMV, EBV)** to host immune dysregulation
* Performed **mediation analyses** connecting microbial signals, host transcriptional modules, and clinical outcomes
* Defined **host endotype dynamics** and their associations with microbial expansion
* Developed **integrated prognostic models** combining host and microbial features

## Repository structure

The repository is organized into modular analysis steps:

```
.
├── 01_overview/        # Cohort overview and descriptive statistics
├── 02_transcriptome/   # Host transcriptome analysis
├── 03_microbiome/      # Microbial profiling and host–microbe interaction
├── 04_endotype/        # Host endotypes, trajectory modeling, and prediction
├── Inputs/             # Processed input data (metadata, counts, modules)
├── Outputs/            # Figures and supplementary data
└── README.md
```

### 01_overview

Basic cohort characterization and summary statistics.

* `overview.R` – cohort description and visualization
* `table1.R` – clinical summary table
* `write_data.R` – data preprocessing / export

---

### 02_transcriptome

Host transcriptional analysis.

* `deg.R` – differential expression analysis
* `deg_dynamics.R` – longitudinal DEG analysis
* `immune.R` – immune marker and pathway analysis

Outputs include:

* DEG results
* Dynamic GSEA
* Immune marker summaries

---

### 03_microbiome

Circulating microbial analysis and host–microbe interactions.

Key analyses:

* Microbial prevalence and abundance
* Mortality association (GLM / MaAsLin3)
* Correlation with host transcriptional modules
* Mediation analysis (microbe → host → outcome)
* Survival analysis

Representative scripts:

* `microbe_mortality_association*.R`
* `microbe_trans_cor.R`
* `mediation.R`, `mediation_d.R`
* `survival analysis.R`
* `Spearman.R`

---

### 04_endotype

Host transcriptional endotypes (CTS/SRS), trajectory dynamics, and predictive modeling.

Key components:

* Endotype definition (CTS, SRS)
* Dynamic transition analysis
* Microbe–endotype interactions
* Mortality prediction models

Modeling approaches:

* GLM / GLMNET
* Random Forest
* SVM
* Gradient Boosting

Representative scripts:

* `CTS.R`, `CTS genes.R`
* `CTSdyn microbe.R`
* `auc_*` (model evaluation and comparison)
* `forest.R`

---

## Analysis workflow

The analysis follows a structured pipeline:

1. **Cohort overview**
   Clinical summary and baseline characteristics

2. **Host transcriptome**

   * DEG analysis
   * Pathway enrichment
   * Longitudinal trajectory inference

3. **Microbiome profiling**

   * Microbial detection and prevalence
   * Mortality association
   * Host–microbe correlation

4. **Integration & mediation**

   * Microbe → host module → outcome
   * Time-lagged associations

5. **Endotype modeling**

   * CTS/SRS classification
   * Dynamic transitions

6. **Prediction modeling**

   * Multi-modal models
   * Longitudinal performance evaluation

---

## Key features

* **Longitudinal design** across multiple ICU timepoints
* **Multi-omics integration** (host + microbes)
* **Dynamic trajectory modeling**
* **Mediation framework** linking microbes and host response
* **Clinically relevant prediction models**

---

## Requirements

Main environment:

* R (≥ 4.0)
* Key packages:

  * `tidyverse`
  * `data.table`
  * `ggplot2`
  * `survival`
  * `glmnet`
  * `caret`
  * `MaAsLin3`
  * `GSVA`
  * `mediation`

