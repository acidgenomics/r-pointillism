# SingleCellExperiment Example Data
# Using splatter to generate simulated counts
# 2018-08-21

library(devtools)
library(tidyverse)
library(bcbioSingleCell)
library(splatter)
library(Seurat)
library(Matrix)
load_all()



# splatter =====================================================================
# note: these DE params are natural log scale
params <- newSplatParams()
params <- setParam(params, "de.facLoc", 1)
params <- setParam(params, "de.facScale", .25)
params <- setParam(params, "dropout.type", "experiment")
params <- setParam(params, "dropout.mid", 3)
sce <- splatSimulate(params, group.prob = c(.5, .5), method = "groups")
# Modify colData
colData(sce) <- camel(colData(sce), rownames = TRUE, colnames = TRUE)
sce$cell <- NULL
# Add sampleID and sampleName columns
sce$batch <- as.factor(camel(sce$batch))
sce$group <- as.factor(camel(sce$group))
sce$sampleID <- sce$group
sce$sampleName <- sce$sampleID
stopifnot("sampleName" %in% colnames(colData(sce)))
# Just keep raw counts
assays(sce) <- assays(sce)["counts"]
# Prepare rowRanges
gr <- makeGRangesFromEnsembl("Homo sapiens", release = 92)
rowRanges(sce) <- gr[seq_len(nrow(sce))]
# Sanitize to camel case
colData(sce) <- camel(colData(sce))
metadata(sce) <- list()
sce <- metrics(sce, recalculate = TRUE)
sce <- filterCells(sce, minCellsPerGene = 25)



# seurat_small =================================================================
seurat_small <- as(sce, "seurat") %>%
    convertGenesToSymbols() %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(resolution = seq(from = 0.4, to = 1.2, by = 0.4)) %>%
    RunTSNE(check_duplicates = FALSE) %>%
    # Requires Python `umap-learn` package
    RunUMAP() %>%
    SetAllIdent(id = "res.0.4")
stopifnot("sampleName" %in% colnames(colData(seurat_small)))



# sce_small ====================================================================
# Convert rows (geneName) back to Ensembl IDs (geneID)
seurat_sce <- seurat_small %>%
    as("SingleCellExperiment") %>%
    convertSymbolsToGenes()
stopifnot(identical(dimnames(sce), dimnames(seurat_sce)))
stopifnot("ident" %in% colnames(colData(seurat_sce)))
# Ensure that dimensional reduction data is slotted correctly
stopifnot(identical(
    names(reducedDims(seurat_sce)),
    c("PCA", "TSNE", "UMAP")
))
assays(sce) <- assays(seurat_sce)
colData(sce) <- colData(seurat_sce)
reducedDims(sce) <- reducedDims(seurat_sce)
# Remove genes with all zero counts prior to zinbwave
sce <- filterCells(sce)
# Calculate zero weights using zinbwave (ZINB-WaVE negative binomial fit)
sce <- runZinbwave(sce)
sce_small <- sce



# marker data frames ===========================================================
all_markers_small <- FindAllMarkers(seurat_small)
all_markers_small <- sanitizeSeuratMarkers(
    data = all_markers_small,
    rowRanges = rowRanges(seurat_small)
)
known_markers_small <- tibble(
    cellType = c("cell_type_1", "cell_type_2"),
    geneID = pull(all_markers_small, "geneID")[seq_len(2L)]
) %>%
    group_by(cellType)
known_markers_detected_small <- knownMarkersDetected(
    all = all_markers_small,
    known = known_markers_small
)



# Save =========================================================================
use_data(
    sce_small,
    seurat_small,
    all_markers_small,
    known_markers_small,
    known_markers_detected_small,
    compress = "xz",
    overwrite = TRUE
)
