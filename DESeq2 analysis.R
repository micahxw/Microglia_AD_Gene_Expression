#DESeq2 analysis.R
#This file takes in GSE146639_metadata and the path to the GSE146639 bulk 
#RNA seq files, and performs DESeq2 analysis and GSEA

library(readxl)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

#replace with your project path to set base path for project
base_path <- "path/to/project/folder" 
metadata_path <- file.path(base_path, "GSE146639_metadata.xlsx")
bulk_path <- file.path(base_path, "bulk")
output_path <- base_path #change output path if needed

#Load metadata
GSE146639_metadata<- read_excel(metadata_path)

#Ensure donor_group has levels AD and CTR
GSE146639_metadata$donor_group <- factor(GSE146639_metadata$donor_group, levels = c("CTR", "AD"))

#Load all bulk RNA-seq count files 
bulk_files<- list.files(bulk_path, full.names = TRUE)

#process count data for DESeq2 analysis, filter for common genes
count_list <- lapply(bulk_files, function(file){
  read.table(file, header = FALSE, row.names = 1)})
gene_lists <- lapply(count_list, rownames)
common_genes <- Reduce(intersect, gene_lists)
filtered_count_list <- lapply(count_list, function(x) x[common_genes, , drop = FALSE])

#create expression df for differential expression analysis
count_matrix <- do.call(cbind, filtered_count_list)
colnames(count_matrix) <- GSE146639_metadata$Sample_id
count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "numeric" 
GSE146639_metadata <- GSE146639_metadata[match(colnames(count_matrix), GSE146639_metadata$Sample_id), ]
rownames(GSE146639_metadata) <- GSE146639_metadata$Sample_id

#Create DESeq2 dataset- adjust for sequencing batch + brain region
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = GSE146639_metadata,
  design = ~ seq_pool + brain_region + donor_group)

#filtering out low-count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

#Run DESeq2
dds <- DESeq(dds)
res <- results(dds, contrast = c("donor_group", "AD", "CTR"))

#Filter and rank by adjusted p-val and effect size
res <- na.omit(res)
ranked_res <- res[order(res$padj), ]
ranked_res <- ranked_res[abs(ranked_res$log2FoldChange) > 0.5, ]
deg_features <- rownames(ranked_res)
# Save to CSV
write.csv(deg_features, file.path(output_path, "deg_features.csv"), row.names = FALSE)

#GSEA for pathway enrichment analysis
#prep DEG list for GSEA
ranked_genes <- res[order(res$log2FoldChange, decreasing = TRUE), ]
gene_ranks <- ranked_genes$log2FoldChange
names(gene_ranks) <- rownames(ranked_genes)

#Run GSEA
gsea_results <- gseGO(
  geneList = gene_ranks,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500
)
#save GSEA results-optional
#write.csv(gsea_results@result, file.path(output_path, "GSEA_results.csv"))

#create dot plot of GSEA results
ggplot(gsea_results, aes(x = reorder(Description, NES), y = NES, size = setSize, color = p.adjust)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Pathway", y = "Normalized Enrichment Score (NES)", 
       title = "Enriched Pathways in AD vs. CTR Microglia", color = "Adjusted p-value") +
  theme_minimal()
