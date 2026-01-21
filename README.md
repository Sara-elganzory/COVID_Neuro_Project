# Unraveling the Transcriptomic Overlap Between SARS-CoV-2 and Neurodegeneration

**Student Name:** Sara Sayed Elganzory
**Student ID:** 221000509
**Course:** Advanced Programming & Data Visualization

![Status](https://img.shields.io/badge/Status-Completed-success)
![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![Focus](https://img.shields.io/badge/Bioinformatics-Transcriptomics-green)

## 📌 Project Overview
This project investigates the molecular mechanisms linking **COVID-19** to the acceleration of **Alzheimer’s Disease (AD)** and **Parkinson’s Disease (PD)**. 

Unlike traditional analyses using web tools (e.g., GEO2R), this project utilizes a **custom Python bioinformatics pipeline** to process raw RNA-Seq data from the Gene Expression Omnibus (GEO). The goal was to identify shared inflammatory and neurodegenerative pathways and to test the "Prion Hypothesis" of SARS-CoV-2.

## 🧬 Key Findings

### 1. COVID-19 & Parkinson's Disease (Oxidative Stress Link)
* **Shared Genes:** Identified **272** genes differentially expressed in both conditions.
* **Key Discovery:** The presence of **MPO (Myeloperoxidase)** and **LPO** in the overlap suggests that COVID-19 may accelerate Parkinsonian pathology through **oxidative stress** and dopaminergic neuron damage.


### 2. COVID-19 & Alzheimer's Disease (Synaptic Link)
* **Shared Genes:** Identified **320** genes differentially expressed in both conditions.
* **Key Discovery:** The dysregulation of **DOC2A** (Synaptic vesicle fusion) and **AGPAT2** indicates a convergence on **synaptic dysfunction** and lipid metabolism.

### 3. The Prion Hypothesis
* **Method:** Screened shared signatures for amyloidogenic markers (*PRNP, SNCA, MAPT, APP*).
* **Outcome:** **No direct transcriptional overlap** was found.
* **Conclusion:** The link between COVID-19 and neurodegeneration is likely driven by systemic inflammation (cytokine storm) rather than direct prion-like protein seeding at the transcriptional level.

## 📂 Repository Structure

```text
COVID_Neuro_Project/
│
├── analysis.py           # The Master Python script (Data processing, Statistics, Plotting)
├── README.md             # Project documentation (You are here)
│
└── results/              # Output folder containing visualizations
    ├── PCA_COVID19.png           # Quality Control Plots
    ├── Volcano_Alzheimers.png    # Differential Expression Plots
    ├── Heatmap_Parkinsons.png    # Hierarchical Clustering
    ├── Venn_COVID_AD.png         # Intersection Analysis
    └── SHARED_COVID_PD.csv       # List of common genes

🛠️ Methodology & Tech Stack
This analysis replaced standard GUI tools with a programmatic approach for reproducibility:

Data Processing: Pandas (Log2 normalization, data cleaning).

Statistical Analysis: SciPy (Welch’s T-Test for unequal variance).

Dimensionality Reduction: Scikit-Learn (Principal Component Analysis).

Visualization: Matplotlib and Seaborn (Volcano plots, Heatmaps).

Intersection Logic: Matplotlib-Venn (Pairwise set analysis).

🚀 How to Run
1. Clone the repository.

2. Install dependencies:

pip install pandas numpy scipy matplotlib seaborn scikit-learn matplotlib-venn mygene

3. Run the analysis script:

python analysis.py

📚 Data Sources

COVID-19: GEO Accession GSE157103

Alzheimer's: GEO Accession GSE159699

Parkinson's: GEO Accession GSE68719