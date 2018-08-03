# seurat_small =================================================================
# Let's handoff to seurat to perform dimensionality reduction and clustering,
# then slot the DR data in our bcbioRNASeq object
seurat_small <- as(sce, "seurat") %>%
    convertGenesToSymbols() %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(resolution = seq(from = 0.4, to = 1.2, by = 0.4)) %>%
    RunTSNE() %>%
    # Requires Python `umap-learn` package
    RunUMAP() %>%
    SetAllIdent(id = "res.0.4")



# all_markers_small ============================================================
all_markers_small <- seurat_small %>%
    FindAllMarkers() %>%
    # Sanitize, including more robust gene annotation information
    sanitizeSeuratMarkers(rowRanges = rowRanges(seurat_small))



# known_markers_small ==========================================================
known_markers_small <- knownMarkersDetected(
    all = all_markers_small,
    known = cell_type_markers[["homoSapiens"]]
)



# cellranger_small =============================================================
# Convert rows (geneName) back to Ensembl IDs (geneID)
seurat_sce <- seurat_small %>%
    as("SingleCellExperiment") %>%
    convertSymbolsToGenes()
# Ensure that dimensional reduction data is slotted correctly
stopifnot(identical(
    names(reducedDims(seurat_sce)),
    c("PCA", "TSNE", "UMAP")
))
stopifnot(identical(dimnames(sce), dimnames(seurat_sce)))
colData(sce) <- colData(seurat_sce)
reducedDims(sce) <- reducedDims(seurat_sce)
cellranger_small <- sce



# Save =========================================================================
use_data(
    cellranger_small,
    seurat_small,
    all_markers_small,
    known_markers_small,
    compress = "xz",
    overwrite = TRUE
)
