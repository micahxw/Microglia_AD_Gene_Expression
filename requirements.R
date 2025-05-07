#This file installs all necessary packages
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "enrichplot"))

# Install CRAN packages
install.packages(c(
  "tidyverse",
  "readxl",
  "writexl",
  "ggplot2",
  "caret",
  "randomForest",
  "e1071",
  "readr",
  "pROC"
))
