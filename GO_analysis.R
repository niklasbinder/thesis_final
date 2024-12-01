BiocManager::install("org.Hs.eg.db", update = FALSE)
BiocManager::install(c(
  "GenomeInfoDbData",
  "AnnotationDbi",
  "GenomeInfoDb",
  "IRanges",
  "S4Vectors"
), update = FALSE)



library(tidyverse)
library(org.Hs.eg.db)

library(clusterProfiler)

de <- read_csv("celltype_differential_expression.csv")
glimpse(de)





macrophage_data <- de %>%
  filter(grepl("Macrophage", cell_type)) %>%
  group_by(cell_type) %>%
  filter(padj   < 0.05 & abs(log2FoldChange) > 1) %>%  
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice_head(n = 50) %>%
  ungroup()


gene_lists <- macrophage_data %>%
  group_by(cell_type) %>%
  summarise(genes = list(gene)) %>%
  deframe()

entrez_lists <- lapply(gene_lists, function(genes) {
  bitr(genes, 
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = org.Hs.eg.db)$ENTREZID
})

compare_go_detailed <- compareCluster(
  geneClusters = entrez_lists,
  OrgDb = org.Hs.eg.db,
  fun = "enrichGO",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 3,    
  maxGSSize = 500  
)



png("monocyte_GO_comparison.png", 
    width = 10, 
    height = 8, 
    units = "in", 
    res = 300,  # This sets it to 300 dpi
    bg = "white")  # White background
dotplot(compare_go_detailed, showCategory = 5, title = "Macrophage GO Biological Process Comparison",font.size = 6)
dev.off()











monocyte_data <- de %>%
  filter(grepl("Monocyte", cell_type)) %>%
  group_by(cell_type) %>%
  filter(padj   < 0.05 & abs(log2FoldChange) > 1) %>%  
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice_head(n = 50) %>%
  ungroup()


gene_lists <- monocyte_data %>%
  group_by(cell_type) %>%
  summarise(genes = list(gene)) %>%
  deframe()

entrez_lists <- lapply(gene_lists, function(genes) {
  bitr(genes, 
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = org.Hs.eg.db)$ENTREZID
})

compare_go_detailed <- compareCluster(
  geneClusters = entrez_lists,
  OrgDb = org.Hs.eg.db,
  fun = "enrichGO",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 3,    
  maxGSSize = 500   
)




  

png("monocyte_GO_comparison.png", 
    width = 10, 
    height = 8, 
    units = "in", 
    res = 300,  # This sets it to 300 dpi
    bg = "white")  # White background

dotplot(compare_go_detailed, showCategory = 5, title = "Monocyte GO Biological Process Comparison",font.size = 7)


dev.off()






