# Microglia Gene Expression in Alzheimer's Disease
# Overview
This project uses bulk RNA-seq data to analyze differential gene expression in microglia in Alzheimer's Disease (AD) vs control (CTR) states. It first identifies top differentially expressed genes (DEGs) in a microglia bulk-RNA dataset, and performs differential expression analysis and GSEA. Then it uses the DEGs as features in a predictive model using general brain bulk-RNA data, in order to explore how microglia contribute to AD.

1. metadata prep.R -> takes in Series Matrices and creates metadata for later analyses (optional, resulting metadata files present in data directory)
2. DESeq2 analysis.R -> takes in GSE146639_metadata and the path to the GSE146639 bulk RNA seq files, and performs DESeq2 analysis and GSEA
3. DEG_modeling.R -> takes in GSE174367 expression data and DEG file from GSE146639 and applies ML models to classify AD vs. CTR
4. requirements.R <- installs all required packages

# Data Preparation:
This project uses **GSE146639** and **GSE174367**, which can be downloaded from the GEO database at 
1. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146639:
   Download the GSE146639_RAW.tar file and (optionally) the Series Matrix file (a preprocessed metadata version is available within the repository). Note that you will have to extract the bulk RNA samples from zipped folders within GSE146639_RAW.tar after download.
   - After selecting and unzipping all the bulk RNA samples, save to folder in your directory named 'bulk'.
   - Make sure to save series matrix as 'GSE146639_series_matrix.xlsx'.
2.  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367:
   Download the GSE174367_bulkRNA_processed.rda.gz file and (optionally) the Series Matrix file (a preprocessed metadata version is available within the repository).
    - Make sure series matrix is named 'GSE174367_bulkRNA_processed.rda.gz'
    - Make sure to save series matrix as 'GSE174367_series_matrix.xlsx'
- After downloading the data, create a folder named 'data' in the root of the repository, and place everything in the folder. Use the folder path in the scripts wherever there's a placeholder path.

# R Package Requirements - Use requirements.R to install
- tidyverse
- readxl
- writexl
- DESeq2
- clusterProfiler
- org.Hs.eg.db
- ggplot2
- caret
- randomForest
- e1071
- readr
- pROC





   
