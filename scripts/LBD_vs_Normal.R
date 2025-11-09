library(DESeq2)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)

# 1. Read input files
counts_file <- read_tsv("counts_clean.tsv")
meta_file   <- read_tsv("samples_metadata.tsv")

# 2. Prepare count matrix (first column = Geneid / Ensembl)
cts <- as.data.frame(counts_file)
rownames(cts) <- cts$Geneid
cts <- cts[, -1]

# 3. Check sample matching between counts and metadata
if(!all(colnames(cts) %in% meta_file$sample)){
  stop("Sample names in counts and metadata do not match. Please check files.")
}

# Reorder metadata to match column order of counts
meta_file <- meta_file[match(colnames(cts), meta_file$sample), ]
rownames(meta_file) <- meta_file$sample

# 4. Create DESeqDataSet (basic model)
dds_basic <- DESeqDataSetFromMatrix(countData = as.matrix(cts),
                                    colData = meta_file,
                                    design = ~ disease_status)

# Filter low-count genes
keep <- rowSums(counts(dds_basic) >= 10) >= 2
dds_basic <- dds_basic[keep, ]

# Run DESeq
dds_basic <- DESeq(dds_basic)

# Extract results: LBD vs Normal
res_basic <- results(dds_basic, contrast = c("disease_status", "LBD", "Normal"))
summary(res_basic)
write.csv(as.data.frame(res_basic), "DESeq2_LBD_vs_Normal_results_basic.csv")

# 5. Sex-adjusted model: design = ~ sex + disease_status
dds_sex <- DESeqDataSetFromMatrix(countData = as.matrix(cts),
                                  colData = meta_file,
                                  design = ~ sex + disease_status)

keep2 <- rowSums(counts(dds_sex) >= 10) >= 2
dds_sex <- dds_sex[keep2, ]

dds_sex <- DESeq(dds_sex)

# Inspect coefficient names and extract disease effect
resultsNames(dds_sex)
res_sexAdj <- results(dds_sex, contrast = c("disease_status", "LBD", "Normal"))
summary(res_sexAdj)
write.csv(as.data.frame(res_sexAdj), "DESeq2_LBD_vs_Normal_sexAdjusted.csv")

# Optional: lfcShrink for more stable log2 fold changes (apeglm)
if(requireNamespace("apeglm", quietly = TRUE)){
  coef_name <- grep("disease_status", resultsNames(dds_sex), value = TRUE)[1]
  res_sexAdj_shr <- lfcShrink(dds_sex, coef = coef_name, type = "apeglm")
  write.csv(as.data.frame(res_sexAdj_shr), "DESeq2_LBD_vs_Normal_sexAdjusted_shrunk.csv")
}

# 6. VST transformation and PCA
vsd <- vst(dds_sex, blind = FALSE)
# PCA plot data
pca_data <- plotPCA(vsd, intgroup = c("disease_status", "sex"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Save PCA plot
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = disease_status, shape = sex)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA: LBD vs Normal (sex-adjusted vst)")
ggsave("pca_LBD_vs_Normal_sexAdjusted.png", pca_plot, width = 7, height = 6)

# 7. Heatmap of top N genes by padj (from sex-adjusted results)
res_df <- as.data.frame(res_sexAdj)
res_df$ensembl <- rownames(res_df)
res_df <- res_df %>% arrange(padj)

topN <- 30
top_genes <- head(rownames(res_sexAdj[order(res_sexAdj$padj), ]), n = topN)
mat <- assay(vsd)[top_genes, ]
mat <- t(scale(t(mat)))  # z-score by gene
annotation_col <- as.data.frame(colData(dds_sex)[, c("disease_status", "sex")])

png("heatmap_top30_LBD_vs_Normal_sexAdjusted.png", width = 1200, height = 1000, res = 150)
pheatmap(mat, annotation_col = annotation_col, show_rownames = TRUE, cluster_cols = TRUE, main = "Top 30 DEGs: LBD vs Normal (sex-adjusted)")
dev.off()

# 8. Volcano plot (basic ggplot)
res_plot_df <- res_df %>%
  mutate(symbol = mapIds(org.Hs.eg.db, keys = ensembl, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")) %>%
  mutate(padj = ifelse(is.na(padj), 1, padj), negLogP = -log10(padj))

res_plot_df <- res_plot_df %>%
  mutate(status = case_when(
    padj < 0.1 & log2FoldChange > 1 ~ "Up",
    padj < 0.1 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "NS"
  ))

volcano <- ggplot(res_plot_df, aes(x = log2FoldChange, y = negLogP, color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Down" = "blue", "NS" = "grey70", "Up" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  theme_minimal() +
  labs(x = "log2 Fold Change (LBD / Normal)", y = "-log10 adjusted p-value", title = "Volcano: LBD vs Normal (sex-adjusted)")

ggsave("Volcano_plot.png", volcano, width = 8, height = 6)

# 9. Functional enrichment (ORA) using clusterProfiler (GO BP)
# Map significant genes (padj < 0.1) to ENTREZ
sig_tbl <- res_df %>% filter(!is.na(padj), padj < 0.1)
sig_map <- bitr(sig_tbl$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Build universe
univ_tbl <- res_df %>% filter(!is.na(symbol))
univ_map <- bitr(unique(univ_tbl$symbol), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene = unique(sig_map$ENTREZID),
                universe = unique(univ_map$ENTREZID),
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.25,
                qvalueCutoff = 0.25,
                readable = TRUE)

ego_df <- as.data.frame(ego)
write.csv(ego_df, "GO_BP_enrich_results_relaxed.csv", row.names = FALSE)

# 10. KEGG enrichment (if available)
ekegg <- enrichKEGG(gene = unique(sig_map$ENTREZID), organism = "hsa", pvalueCutoff = 0.25, qvalueCutoff = 0.25)
ekegg_df <- as.data.frame(ekegg)
write.csv(ekegg_df, "KEGG_enrich_results_relaxed.csv", row.names = FALSE)

