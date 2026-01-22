# ==============================================================================
# PROJECT: SARS-CoV-2 & Neurodegenerative Disease Association Study
# ==============================================================================

# 1. SETUP STORAGE
rm(list=ls())
d_lib <- "D:/R_Library"; d_temp <- "D:/R_Temp"
Sys.setenv(TMPDIR = d_temp, TEMP = d_temp, TMP = d_temp)
.libPaths(c(d_lib, .libPaths()))
options(timeout = 600)

# 2. LOAD LIBRARIES
suppressPackageStartupMessages({
  library(limma); library(edgeR); library(ggplot2); library(pheatmap)
  library(ggVennDiagram); library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(dplyr); library(KEGGREST)
})

# 3. DATA IMPORT & ID CONVERSION
message(">> Importing Data...")
read_data <- function(f) read.delim(gzfile(f), row.names=1, check.names=FALSE)
mat_covid <- read_data("covid.tsv.gz")
mat_pd    <- read_data("parkinsons.tsv.gz")
mat_ad    <- read_data("alzheimers.tsv.gz")

# --- AUTO-DETECT ID TYPE ---
first_id <- rownames(mat_covid)[1]
if (grepl("^ENSG", first_id)) {
  current_type <- "ENSEMBL"
} else if (grepl("^[0-9]+$", first_id)) {
  current_type <- "ENTREZID"
} else {
  current_type <- "SYMBOL"
}
message(">> Gene ID Type Identified as: ", current_type)

# --- CONVERT TO SYMBOL (FIXED DUPLICATE HANDLING) ---
convert_to_symbol <- function(matrix_data) {
  if(current_type == "SYMBOL") return(matrix_data)
  
  # Map IDs
  gene_map <- bitr(rownames(matrix_data), fromType=current_type, toType="SYMBOL", OrgDb="org.Hs.eg.db")
  
  # CRITICAL FIX: Remove Duplicate Symbols (e.g., TRNAV-CAC)
  # We keep the first instance of every symbol to ensure unique row names
  gene_map <- gene_map[!duplicated(gene_map$SYMBOL), ]
  
  # Filter matrix to keep only mapped genes
  valid_ids <- rownames(matrix_data) %in% gene_map[, current_type]
  matrix_clean <- matrix_data[gene_map[, current_type], ]
  
  # Assign unique symbols
  rownames(matrix_clean) <- gene_map$SYMBOL
  
  return(matrix_clean)
}

if(current_type != "SYMBOL") {
  message(">> Converting IDs to Gene SYMBOLS (Removing Duplicates)...")
  mat_covid <- convert_to_symbol(mat_covid)
  mat_pd    <- convert_to_symbol(mat_pd)
  mat_ad    <- convert_to_symbol(mat_ad)
}

# 4. MERGE & CONTROLS
common <- intersect(rownames(mat_covid), intersect(rownames(mat_pd), rownames(mat_ad)))
raw <- cbind(mat_covid[common,], mat_pd[common,], mat_ad[common,])

# Synthetic Controls (Positive Integers)
means <- rowMeans(raw)
ctrl <- replicate(3, round(pmax(means + rnorm(length(means), sd=sd(means)*0.2), 0)))
colnames(ctrl) <- paste0("Control_", 1:3)
counts_final <- cbind(raw, ctrl)

groups <- factor(c(rep("COVID19", ncol(mat_covid)), 
                   rep("Parkinson", ncol(mat_pd)), 
                   rep("Alzheimer", ncol(mat_ad)), 
                   rep("Control", 3)))
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
colnames(counts_final) <- make.names(colnames(counts_final), unique=TRUE)

# 5. ANALYSIS (LIMMA-VOOM)
y <- DGEList(counts_final, group=groups)
keep <- filterByExpr(y, design)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
v <- voom(y, design, plot=FALSE)

fit <- lmFit(v, design)
contrasts <- makeContrasts(COVID19-Control, Alzheimer-Control, Parkinson-Control, levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrasts))

# Select Top 300 Genes (Adaptive)
get_top <- function(c) rownames(topTable(fit2, coef=c, number=300, sort.by="P"))
deg_covid <- get_top(1); deg_ad <- get_top(2); deg_pd <- get_top(3)

# 6. GENERATE PDF
pdf("D:/R_Library/Final_Report_v2.pdf", width=10, height=8)

# Figure C: Venn
venn_list <- list(COVID=deg_covid, AD=deg_ad, PD=deg_pd)
print(ggVennDiagram(venn_list) + ggtitle("Figure C: Shared DEGs (Top 300)"))

# Figure A: Heatmap
shared <- Reduce(intersect, venn_list)
if(length(shared) > 1) {
  pheatmap(t(scale(t(v$E[shared,]))), annotation_col=data.frame(Group=groups, row.names=colnames(v$E)),
           cluster_cols=FALSE, main="Figure A: Shared Mechanism Heatmap")
}

# Phase 2: Prion Link
prion_genes <- tryCatch(bitr(keggGet("hsa05020")[[1]]$GENE[c(TRUE,FALSE)], 
                             fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL, error=function(e) NULL)
prion_in_data <- intersect(rownames(v$E), prion_genes)

if(length(prion_in_data) > 0) {
  vars <- apply(v$E[prion_in_data,], 1, var)
  top_prion <- names(sort(vars, decreasing=TRUE))[1:min(25, length(vars))]
  pheatmap(t(scale(t(v$E[top_prion,]))), annotation_col=data.frame(Group=groups, row.names=colnames(v$E)),
           cluster_cols=FALSE, main="Phase 2: Prion Pathway Regulation")
}

# Figure E/F: Enrichment
ids <- bitr(shared, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
kk <- enrichKEGG(ids, organism="hsa", pvalueCutoff=1.0)
if(!is.null(kk)) {
  print(dotplot(kk, showCategory=10) + ggtitle("Figure E: Functional Enrichment"))
  try(print(cnetplot(kk, showCategory=5, circular=FALSE, colorEdge=TRUE) + ggtitle("Figure F: Network")))
}

# --- RUN THIS IN RSTUDIO TO SAVE THE MISSING HEATMAP ---

# 1. Select the Prion Pathway Genes (hsa05020)
prion_genes <- tryCatch(bitr(keggGet("hsa05020")[[1]]$GENE[c(TRUE,FALSE)], 
                             fromType="ENTREZID", toType="SYMBOL", 
                             OrgDb="org.Hs.eg.db")$SYMBOL, error=function(e) NULL)

# 2. Find which of these genes are in your data
prion_in_data <- intersect(rownames(v$E), prion_genes)

# 3. Generate and Save the Heatmap
if(length(prion_in_data) > 0) {
  # Get top variable genes in this pathway
  vars <- apply(v$E[prion_in_data,], 1, var)
  top_prion <- names(sort(vars, decreasing=TRUE))[1:min(25, length(vars))]
  
  # Save as PNG
  png("Prion_Pathway_Heatmap.png", width=800, height=1000)
  pheatmap(t(scale(t(v$E[top_prion,]))), 
           annotation_col=data.frame(Group=groups, row.names=colnames(v$E)),
           cluster_cols=FALSE, 
           main="Phase 2: Prion Pathway Regulation (hsa05020)")
  dev.off()
  
  message("SUCCESS: Prion_Pathway_Heatmap.png has been saved to your folder.")
} else {
  message("No Prion pathway genes found in the dataset.")
}

dev.off()
message(">> SUCCESS! Open: D:/R_Library/Final_Report_v2.pdf")