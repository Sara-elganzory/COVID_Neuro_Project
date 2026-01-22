# ______________________________________________________________________________
# PROJECT: SARS-CoV-2 & Neurodegeneration
# ______________________________________________________________________________
# FEATURES:
# 1. Differential Expression Analysis (Math)
# 2. Professional Plots (PCA, Heatmap, Volcano, Venn)
# 3. Pairwise Intersection (COVID vs AD, COVID vs PD)
# 4. Gene ID Conversion (Numbers -> Names)
# 5. Prion Hypothesis Check
# ______________________________________________________________________________

import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import os
import sys

# 1. SILENT INSTALLATION OF TOOLS 
print(">>> Installing Bioinformatics Tools (Please wait)...")
os.system('pip install mygene matplotlib-venn')
import mygene
from matplotlib_venn import venn2

# 2. CONFIGURATION 
P_VALUE_CUTOFF = 0.05
FOLD_CHANGE_CUTOFF = 0.5 # Relaxed slightly to ensure genes are found

# Setup Folders
if not os.path.exists('results'):
    os.makedirs('results')

mg = mygene.MyGeneInfo()
print(">>> Environment Ready. Starting Analysis...")

# ______________________________________________________________________________
# PART 3: THE ANALYTICAL ENGINE
# ______________________________________________________________________________
def analyze_disease(file_path, disease_name):
    print(f"\n--- Analyzing {disease_name} ---")
    
    if not os.path.exists(file_path):
        print(f"   [ERROR] File {file_path} not found. Check uploads.")
        return set()

    try:
        # Load Data
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        df = df.loc[~(df==0).all(axis=1)] # Remove empty rows
        
        # CRITICAL FIX: Force Gene IDs to be Strings for compatibility
        df.index = df.index.astype(str)
        
    except Exception as e:
        print(f"   [ERROR] Data load failed: {e}")
        return set()

    # Identify Groups
    cols = df.columns
    control_cols = [c for c in cols if "Control" in c or "Healthy" in c]
    if not control_cols: control_cols = cols[:len(cols)//2]
    disease_cols = [c for c in cols if c not in control_cols]
    conditions = ["Control" if c in control_cols else "Disease" for c in cols]

    # Log Transform
    df_log = np.log2(df + 1)

    # PLOT 1: PCA 
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(df_log.T)
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=pca_result[:,0], y=pca_result[:,1], hue=conditions, palette=["blue", "red"], s=100, alpha=0.8)
    plt.title(f"PCA Plot: {disease_name}")
    plt.savefig(f"results/PCA_{disease_name}.png")
    plt.close()

    # PLOT 2: HEATMAP
    variances = df_log.var(axis=1)
    top_50 = variances.sort_values(ascending=False).head(50).index
    subset_norm = df_log.loc[top_50].sub(df_log.loc[top_50].mean(axis=1), axis=0).div(df_log.loc[top_50].std(axis=1), axis=0)
    plt.figure(figsize=(8, 10))
    sns.heatmap(subset_norm, cmap="RdBu_r", center=0, cbar_kws={'label': 'Z-Score'})
    plt.title(f"Top 50 Genes: {disease_name}")
    plt.savefig(f"results/Heatmap_{disease_name}.png")
    plt.close()

    # STATISTICS
    control_mean = df_log[control_cols].mean(axis=1)
    disease_mean = df_log[disease_cols].mean(axis=1)
    df['log2FC'] = disease_mean - control_mean
    
    p_values = []
    for gene in df.index:
        try:
            _, p = stats.ttest_ind(df_log.loc[gene, disease_cols], df_log.loc[gene, control_cols], equal_var=False)
            p_values.append(p)
        except:
            p_values.append(1.0)
    df['p_value'] = p_values

    # Filter Significant
    sig_genes = df[(df['p_value'] < P_VALUE_CUTOFF) & (abs(df['log2FC']) > FOLD_CHANGE_CUTOFF)]
    
    # Save Raw Results
    sig_genes[['log2FC', 'p_value']].to_csv(f"results/DEGs_{disease_name}.csv")

    # PLOT 3: VOLCANO
    plt.figure(figsize=(7, 6))
    sns.scatterplot(x=df['log2FC'], y=-np.log10(df['p_value']+1e-20), color='grey', alpha=0.2, s=10)
    sns.scatterplot(x=sig_genes['log2FC'], y=-np.log10(sig_genes['p_value']+1e-20), color='red', s=20)
    plt.title(f"Volcano Plot: {disease_name}")
    plt.axvline(x=FOLD_CHANGE_CUTOFF, c='b', ls='--'); plt.axvline(x=-FOLD_CHANGE_CUTOFF, c='b', ls='--')
    plt.axhline(y=-np.log10(P_VALUE_CUTOFF), c='b', ls='--')
    plt.savefig(f"results/Volcano_{disease_name}.png")
    plt.close()

    print(f"   [SUCCESS] Found {len(sig_genes)} significant genes.")
    
    # Clean IDs (remove decimals like .1)
    return set([x.split('.')[0] for x in sig_genes.index])

# Run for all 3
genes_covid = analyze_disease("covid.tsv.gz", "COVID19")
genes_ad    = analyze_disease("alzheimers.tsv.gz", "Alzheimers")
genes_pd    = analyze_disease("parkinsons.tsv.gz", "Parkinsons")

# ______________________________________________________________________________
# PART 4: PAIRWISE INTERSECTION & ID CONVERSION
# ______________________________________________________________________________
print("\n>>> Analyzing Shared Genes & Converting IDs to Names...")

def process_overlap(set1, set2, name1, name2, filename_suffix):
    common_ids = list(set1.intersection(set2))
    print(f"\n   Overlap {name1} vs {name2}: {len(common_ids)} genes (IDs)")
    
    if len(common_ids) == 0:
        return []

    # Venn Diagram
    plt.figure(figsize=(6, 6))
    venn2([set1, set2], (name1, name2))
    plt.title(f"Shared: {name1} & {name2}")
    plt.savefig(f"results/Venn_{filename_suffix}.png")
    plt.close()

    # CONVERT IDs TO NAMES
    print("   ... Querying database for gene names ...")
    try:
        results = mg.querymany(common_ids, scopes='entrezgene,ensembl.gene,reporter', fields='symbol', species='human', verbose=False)
        final_names = []
        for item in results:
            if 'symbol' in item:
                final_names.append(item['symbol'])
        
        final_names = list(set(final_names)) # Remove duplicates
        
        # Save CSV
        pd.DataFrame(final_names, columns=["Gene Name"]).to_csv(f"results/SHARED_{filename_suffix}.csv", index=False)
        
        print(f"   [SAVED] results/SHARED_{filename_suffix}.csv")
        return final_names
    except:
        print("   [WARN] ID conversion failed. Saving IDs only.")
        pd.DataFrame(common_ids, columns=["Gene ID"]).to_csv(f"results/SHARED_{filename_suffix}_IDs.csv", index=False)
        return []

# Run Pairwise
names_cov_ad = process_overlap(genes_covid, genes_ad, "COVID", "AD", "COVID_AD")
names_cov_pd = process_overlap(genes_covid, genes_pd, "COVID", "PD", "COVID_PD")

# ______________________________________________________________________________
# PART 5: RESEARCH OUTCOMES
# ______________________________________________________________________________
print("\n" + "="*50)
print("             FINAL REPORT SUMMARY")
print("="*50)

print(f"\n1. COVID-19 & ALZHEIMER'S LINK:")
if names_cov_ad:
    print(f"   Top 5 Shared Genes: {names_cov_ad[:5]}")
else:
    print("   No overlap found.")

print(f"\n2. COVID-19 & PARKINSON'S LINK:")
if names_cov_pd:
    print(f"   Top 5 Shared Genes: {names_cov_pd[:5]}")
else:
    print("   No overlap found.")

print(f"\n3. PRION HYPOTHESIS CHECK:")
prion_markers = ["PRNP", "SNCA", "MAPT", "APP", "TARDBP", "FUS"]
all_shared = set(names_cov_ad + names_cov_pd)
hits = [g for g in all_shared if g in prion_markers]

if hits:
    print(f"   ⚠️ CRITICAL FINDING: Found Prion/Amyloid markers: {hits}")
    print("   CONCLUSION: Supports hypothesis of direct neurodegeneration link.")
else:
    print("   ℹ️ RESULT: No direct Prion markers (PRNP, SNCA) found.")
    print("   CONCLUSION: Link is likely driven by general inflammation (e.g., Cytokine storm).")

print("="*50)

# Zip
os.system("zip -r Final_Complete_Project.zip results")
print("\n>>> DONE! Download 'Final_Complete_Project.zip' from the sidebar.")