# Functional enrichment and pathway analysis using GO, KEGG, and REACTOME

# Packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(circlize)
library(RColorBrewer)

# Protein list
genes_symbol <- c("CHI3L1", "CYCS", "DDAH1", "GDA", "LRRC4B", "NPTX2", "PKM", "SMOC1", "SPON1", "YWHAG", "YWHAZ")

# GO enrichment (BP, CC, MF)
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
  df <- df[sapply(strsplit(df$geneID, "/"), length) >= 2, ]
  df$code <- paste0(category, seq_len(nrow(df)))
  data.frame(Category = category,
             Code = df$code,
             ID = df$ID,
             term = df$Description,
             genes = df$geneID,
             adj_pval = df$p.adjust)
}
terms_go <- do.call(rbind, Map(prep_terms, go_results, ontologies))

# Conversion to Entrez ID
genes_entrez_df <- bitr(genes_symbol, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Hs.eg.db)
genes_entrez <- genes_entrez_df$ENTREZID

# KEGG and Reactome enrichment
ekegg <- enrichKEGG(gene = genes_entrez, organism = "hsa", pvalueCutoff = 0.05)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
ereact <- enrichPathway(gene = genes_entrez, organism = "human",
                        readable = TRUE, pvalueCutoff = 0.05)
prep_pathway <- function(obj, cat) {
  df <- obj@result
  df <- subset(df, p.adjust <= 0.05)
  df <- df[sapply(strsplit(df$geneID, "/"), length) >= 2, ]
  df$code <- paste0(cat, seq_len(nrow(df)))
  data.frame(Category = cat,
             Code = df$code,
             ID = df$ID,
             term = df$Description,
             genes = df$geneID,
             adj_pval = df$p.adjust)
}
terms_kegg   <- prep_pathway(ekegg, "KEGG")
terms_react  <- prep_pathway(ereact, "RP")

# Save results
write.csv(terms_go, "GOenrichmentPPAV11.csv", row.names = FALSE)
write.csv(rbind(terms_kegg, terms_react), "Pathway_enrichment_filtered.csv", row.names = FALSE)

# Fixed color palette per protein
colors_genes <- setNames(c("olivedrab3", "darkorchid4", "blue4",
                            "aquamarine3", "darkkhaki", "cyan3",
                            "plum3", "violet", "darkseagreen2",
                            "deeppink2", "gold1", "cornflowerblue", "chocolate1",
                            "coral2", "maroon4", "darkmagenta", "darkcyan",
                            "darkred", "chartreuse2", "cadetblue2",
                            "azure2", "yellow1", "hotpink2", "brown4", "mediumpurple3",
                            "bisque2", "dodgerblue2", "forestgreen", "yellowgreen",
                            "skyblue", "violetred", "pink4", "seagreen3")[1:length(genes_symbol)],
                          genes_symbol) 

# GO plot
long_go <- do.call(rbind, lapply(1:nrow(terms_go), function(i) {
  genes <- unlist(strsplit(terms_go$genes[i], "/"))
  genes <- genes[genes %in% names(colors_genes)]
  if (length(genes) == 0) return(NULL)
  data.frame(term = terms_go$Code[i], variable = genes)
}))
chordDiagram(long_go,
             annotationTrack = "grid",
             transparency = 0.5,
             grid.col = colors_genes,
             col = colors_genes[long_go$variable],
             preAllocateTracks = list(track.height = 0.05))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  circos.text(mean(get.cell.meta.data("xlim")),
              get.cell.meta.data("ylim")[2] + mm_y(2),
              get.cell.meta.data("sector.index"),
              facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = 0.6)
}, bg.border = NA)
circos.clear()
print (terms_go)

generate_GO_diagram <- function(genes_symbol, file_name = "GOdiagram.pdf") {
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
    df <- df[sapply(strsplit(df$geneID, "/"), length) >= 2, ]
    df$code <- paste0(category, seq_len(nrow(df)))
    data.frame(Category = category,
               Code = df$code,
               ID = df$ID,
               term = df$Description,
               genes = df$geneID,
               adj_pval = df$p.adjust)
  }
  terms_go <- do.call(rbind, Map(prep_terms, go_results, ontologies))
  enriched_genes <- unique(unlist(strsplit(terms_go$genes, "/")))
  colors_genes <- setNames(colorRampPalette(brewer.pal(8, "Set3"))(length(enriched_genes)),
                            enriched_genes)
  long_go <- do.call(rbind, lapply(1:nrow(terms_go), function(i) {
    genes <- unlist(strsplit(terms_go$genes[i], "/"))
    data.frame(term = terms_go$Code[i], variable = genes)
  }))
  long_go <- subset(long_go, variable %in% names(colors_genes))

  if (nrow(long_go) > 0) {
    pdf(file_name, width = 8, height = 8)
    chordDiagram(long_go,
                 annotationTrack = "grid",
                 transparency = 0.5,
                 grid.col = colors_genes,
               col = colors_genes[long_go$variable],
               preAllocateTracks = list(track.height = 0.05))
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    circos.text(mean(get.cell.meta.data("xlim")),
                get.cell.meta.data("ylim")[2] + mm_y(2),
                get.cell.meta.data("sector.index"),
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.6)
  }, bg.border = NA)
  circos.clear()
  dev.off()
} else {
  message("No valid GO terms with at least 2 genes to plot.")
}

return(terms_go)
}

# KEGG and Reactome plot
terms_pathways <- rbind(terms_kegg, terms_react)
long_path <- do.call(rbind, lapply(1:nrow(terms_pathways), function(i) {
  genes <- unlist(strsplit(terms_pathways$genes[i], "/"))
  genes <- genes[genes %in% names(colors_genes)]
  if (length(genes) == 0) return(NULL)
  data.frame(term = terms_pathways$Code[i], variable = genes)
}))

chordDiagram(long_path,
             annotationTrack = "grid",
             transparency = 0.5,
             grid.col = colors_genes,
             col = colors_genes[long_path$variable],
             preAllocateTracks = list(track.height = 0.05))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  circos.text(mean(get.cell.meta.data("xlim")),
              get.cell.meta.data("ylim")[2] + mm_y(2),
              get.cell.meta.data("sector.index"),
              facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = 0.6)
}, bg.border = NA)
circos.clear()
print(terms_pathways)

enriched_genes <- unique(unlist(strsplit(terms_pathways$genes, "/")))
colors_genes <- setNames(colorRampPalette(c("olivedrab3", "darkorchid4", "blue4",
                                            "aquamarine3", "maroon4", "cyan3",
                                            "plum3", "violet", "darkseagreen2",
                                            "deeppink2", "gold1", "tomato", "tan"))(length(enriched_genes)),
                         enriched_genes)

long_path <- do.call(rbind, lapply(1:nrow(terms_pathways), function(i) {
  genes <- unlist(strsplit(terms_pathways$genes[i], "/"))
  data.frame(term = terms_pathways$Code[i], variable = genes)
}))
long_path <- subset(long_path, variable %in% names(colors_genes))

if (nrow(long_path) > 0) {
  chordDiagram(long_path,
               annotationTrack = "grid",
               transparency = 0.5,
               grid.col = colors_genes,
               col = colors_genes[long_path$variable],
               preAllocateTracks = list(track.height = 0.05))
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    circos.text(mean(get.cell.meta.data("xlim")),
                get.cell.meta.data("ylim")[2] + mm_y(2),
                get.cell.meta.data("sector.index"),
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.6)
  }, bg.border = NA)
  circos.clear()
} else {
  message("No valid genes to plot.")
}

print(terms_pathways)

# Save results
write.csv(rbind(terms_go[, c("Code", "term")],
                terms_kegg[, c("Code", "term")],
                terms_react[, c("Code", "term")]),
          "Term_Code_Reference.csv", row.names = FALSE)
