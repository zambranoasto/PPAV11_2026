# Functional enrichment and pathway analysis of PPAV11 using GO, KEGG, and REACTOME databases

# Import required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(circlize)
library(RColorBrewer)

# Define protein list
genes_symbol <- c("CHI3L1", "CYCS", "DDAH1", "GDA", "LRRC4B", "NPTX2", "PKM", "SMOC1", "SPON1", "YWHAG", "YWHAZ")

# Run GO enrichment analyses (BP, CC, MF)
ontologies <- c("BP", "CC", "MF")
go_results <- lapply(ontologies, function(ont) {
  enrichGO(gene = genes_symbol, OrgDb = org.Hs.eg.db,
           keyType = "SYMBOL", ont = ont,
           pAdjustMethod = "BH", pvalueCutoff = 0.05)
})
names(go_results) <- ontologies

prep_terms <- function(res, category) {
  df <- res@result
  df <- subset(df, p.adjust <= 0.05)
  df <- df[sapply(strsplit(as.character(df$geneID), "/"), length) >= 2, ]
  if (nrow(df) == 0) return(NULL) 
  df$code <- paste0(category, seq_len(nrow(df)))
  data.frame(Category = category,
             Code = df$code,
             ID = df$ID,
             term = df$Description,
             genes = df$geneID,
             adj_pval = df$p.adjust)
}
terms_go <- do.call(rbind, Map(prep_terms, go_results, ontologies))

# Run KEGG and REACTOME enrichment analyses
genes_entrez_df <- bitr(genes_symbol, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Hs.eg.db)
genes_entrez <- genes_entrez_df$ENTREZID
ekegg <- enrichKEGG(gene = genes_entrez, organism = "hsa", pvalueCutoff = 0.05)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
ereact <- enrichPathway(gene = genes_entrez, organism = "human",
                        readable = TRUE, pvalueCutoff = 0.05)
terms_kegg  <- prep_terms(ekegg, "KEGG")
terms_react <- prep_terms(ereact, "RP")
terms_pathways <- rbind(terms_kegg, terms_react)

# Save results
write.csv(terms_go, "Output_GO.csv", row.names = FALSE)
write.csv(rbind(terms_kegg, terms_react), "Output_kegg_reactome.csv", row.names = FALSE)

# Define color pallete
all_enriched_genes <- unique(unlist(strsplit(as.character(c(terms_go$genes, terms_pathways$genes)), "/")))
colors_genes <- setNames(colorRampPalette(brewer.pal(8, "Set3"))(length(all_enriched_genes)),
                         all_enriched_genes)

# Generate GO plot
if (!is.null(terms_go) && nrow(terms_go) > 0) {
  # Force a new plot window in R
  dev.new() 
  long_go <- do.call(rbind, lapply(1:nrow(terms_go), function(i) {
    genes <- unlist(strsplit(as.character(terms_go$genes[i]), "/"))
    data.frame(term = terms_go$Code[i], variable = genes)
  }))
  chordDiagram(long_go,
               annotationTrack = "grid",
               transparency = 0.5,
               grid.col = colors_genes,
               col = colors_genes[as.character(long_go$variable)],
               preAllocateTracks = list(track.height = 0.05))
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    circos.text(mean(get.cell.meta.data("xlim")),
                get.cell.meta.data("ylim")[2] + mm_y(2),
                get.cell.meta.data("sector.index"),
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.6)
  }, bg.border = NA)
  title("Gene Ontology Analysis")
  circos.clear()
}

# Generate KEGG and REACTOME plot
if (!is.null(terms_pathways) && nrow(terms_pathways) > 0) {
  # Force a second new plot window in R
  dev.new() 
  long_path <- do.call(rbind, lapply(1:nrow(terms_pathways), function(i) {
    genes <- unlist(strsplit(as.character(terms_pathways$genes[i]), "/"))
    data.frame(term = terms_pathways$Code[i], variable = genes)
  }))
  chordDiagram(long_path,
               annotationTrack = "grid",
               transparency = 0.5,
               grid.col = colors_genes,
               col = colors_genes[as.character(long_path$variable)],
               preAllocateTracks = list(track.height = 0.05))
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    circos.text(mean(get.cell.meta.data("xlim")),
                get.cell.meta.data("ylim")[2] + mm_y(2),
                get.cell.meta.data("sector.index"),
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.6)
  }, bg.border = NA)
  
  title("KEGG & Reactome Analysis")
  circos.clear()
}
