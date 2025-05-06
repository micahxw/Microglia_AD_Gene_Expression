library(DESeq2)
library(readxl)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

#replace with your base path here
base_path <- "path/to/project/folder"
metadata_path <- file.path(base_path, "GSE146639_metadata.xlsx")
bulk_path <- file.path(base_path, "bulk")
output_path <- base_path

#Load metadata
GSE146639_metadata<- read_excel(metadata_path)

#Ensure donor_group has levels AD and CTR
GSE146639_metadata$donor_group <- factor(GSE146639_metadata$donor_group, levels = c("CTR", "AD"))

#Load all bulk RNA-seq count files 
bulk_files<- list.files(bulk_path, full.names = TRUE)

#process for DESeq2 analysis
count_list <- lapply(bulk_files, function(file){
  read.table(file, header = FALSE, row.names = 1)})
gene_lists <- lapply(count_list, rownames)
common_genes <- Reduce(intersect, gene_lists)
filtered_count_list <- lapply(count_list, function(x) x[common_genes, , drop = FALSE])

#create expression df for differential expression analysis
count_matrix <- do.call(cbind, filtered_count_list)
colnames(count_matrix) <- GSE146639_metadata$Sample_id

#Save GSE146639 expression matrix for downstream modeling
write.csv(count_matrix, file.path(output_path, "146639_expr_data.csv"))

#continue prepping count matrix and metadata for DESeq2
count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "numeric" 
GSE146639_metadata <- GSE146639_metadata[match(colnames(count_matrix), GSE146639_metadata$Sample_id), ]
rownames(GSE146639_metadata) <- GSE146639_metadata$Sample_id

#Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = GSE146639_metadata,
  design = ~ donor_group)

dds <- DESeq(dds)
res <- results(dds)

#Select top DEGs for ML features
ranked_res <- res[order(res$pvalue), ]
deg_features <- rownames(ranked_res)[1:250]
write.csv(deg_features, file.path(output_path, "deg_features.csv"), row.names = FALSE)


#Save variance-stabilized data for potential downstream analysis
vsd <- vst(dds, blind = FALSE)
norm_counts <- assay(vsd)
write.csv(norm_counts, file.path(output_path, "norm_counts.csv"))
write.csv(rownames(norm_counts), file.path(output_path, "gene_names.csv"), row.names = FALSE)

#GSEA for pathway enrichment analysis
res <- na.omit(res)
ranked_genes <- res[order(res$log2FoldChange, decreasing = TRUE), ]
gene_ranks <- ranked_genes$log2FoldChange
names(gene_ranks) <- rownames(ranked_genes)

gsea_results <- gseGO(
  geneList = gene_ranks,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500
)
#save GSEA results
write.csv(gsea_results@result, file.path(output_path, "GSEA_results.csv"))

#create dot plot of GSEA results
ggplot(gsea_results, aes(x = reorder(Description, NES), y = NES, size = setSize, color = p.adjust)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Pathway", y = "Normalized Enrichment Score (NES)", 
       title = "Dot Plot for GSEA Results", color = "Adjusted p-value") +
  theme_minimal()