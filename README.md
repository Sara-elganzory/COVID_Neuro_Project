# 🧬 Unraveling the Transcriptomic Overlap Between SARS-CoV-2 and Neurodegeneration

**Student Name:** Sara Sayed Elganzory
**Student ID:** 221000509
**Course:** Advanced Programming & Data Visualization

![Status](https://img.shields.io/badge/Status-Completed-success)
![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![R](https://img.shields.io/badge/R-Bioconductor-blue)
![Focus](https://img.shields.io/badge/Bioinformatics-Transcriptomics-green)

## 📌 Project Overview

This project investigates the molecular links between **COVID-19 (SARS-CoV-2)** and neurodegenerative conditions (**Alzheimer's** and **Parkinson's Disease**). Using a **Dual-Pipeline Approach (Python & R)**, we identified shared transcriptomic signatures that suggest COVID-19 may accelerate neurodegeneration through indirect environmental stress rather than direct prion protein activation.

## 🎥 Project Demo
**Watch the full video walkthrough of the analysis and findings:**

[![Watch the video](https://img.shields.io/badge/Watch_Video-Click_Here-red?style=for-the-badge&logo=youtube)]([https://youtu.be/_kP7vMhSVms])
> **Note:** The video covers the workflow, key biological findings (MPO/DOC2A), and the explanation of the Prion Hypothesis results.

## 🧪 Objective

To investigate the molecular mechanisms linking COVID-19 to the acceleration of Alzheimer’s (AD) and Parkinson’s (PD), and to specifically test the **"Prion Hypothesis"** (whether SARS-CoV-2 induces prion-like misfolding) using multi-omics validation.

## 🔬 Methodology: A Dual-Pipeline Approach

To ensure reproducibility and statistical robustness, this study utilized two independent analytical pipelines:

### Phase 1: Python (Exploratory Analysis)
* **Goal:** Data cleaning, dimensionality reduction (PCA), and candidate discovery.
* **Tools:** Pandas, SciPy, Scikit-Learn, Matplotlib.
* **Method:** Custom scripts replacing standard GUI tools (GEO2R) for greater control over normalization and intersection logic.

### Phase 2: R/Bioconductor (Statistical Validation)
* **Goal:** Validation of candidates and advanced pathway enrichment.
* **Tools:** limma, edgeR, clusterProfiler, KEGGREST.
* **Method:** Voom normalization for mean-variance trend correction and Empirical Bayes statistics for differential expression.

## 🧬 Key Biological Insights

### 1. COVID-19 & Parkinson's Disease (The Oxidative Stress Link)
* **Shared Genes:** Identified **272** genes differentially expressed in both conditions.
* **Key Discovery:** The upregulation of **MPO (Myeloperoxidase)** and **LPO** in the overlap.
* **Implication:** COVID-19 triggers the same **oxidative stress** pathways that drive dopaminergic neuron death in Parkinson’s pathology.

### 2. COVID-19 & Alzheimer's Disease (The Synaptic Link)
* **Shared Genes:** Identified **320** genes differentially expressed in both conditions.
* **Key Discovery:** The dysregulation of **DOC2A** (Synaptic vesicle fusion) and **AGPAT2**.
* **Implication:** A convergence on **synaptic dysfunction** and lipid metabolism, potentially explaining the cognitive "brain fog" associated with Long-COVID.

### 3. The Prion Hypothesis (The "Two-Hit" Model)
* **Initial Screen (Python):** No direct transcriptional overlap was found for *PRNP*, *SNCA*, or *MAPT*.
* **Pathway Validation (R):** Targeted enrichment revealed systemic dysregulation of the **Prion Diseases Pathway (hsa05020)**.
* **Conclusion:** The link is **indirect**. SARS-CoV-2 does not directly transcribe prion proteins but induces a cellular stress environment (proteostatic disruption) that mimics the conditions required for prion propagation.

## 📂 Repository Structure

```text
COVID_Neuro_Project/
│
├── README.md                  # Project documentation
│
├── Code/
│   ├── Phase1_Python/         # Exploratory Analysis
│   │   └── analysis_pipeline.py
│   │
│   └── Phase2_R_Validation/   # Statistical Validation
│       └── Final_Project_Code.R
│
├── Results/                   # Output Visualizations
│   ├── PCA_COVID19.png        # Quality Control Plots
│   ├── Venn_COVID_AD.png      # Intersection Analysis
│   ├── Prion_Heatmap.png      # Pathway Dysregulation (R)
│   └── Shared_Genes.csv       # Final List of Common Genes
│
└── Report/
    └── Final_Project_Report.pdf




