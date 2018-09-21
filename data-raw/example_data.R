# SingleCellExperiment Example Data
# 2018-09-21

library(reticulate)
library(bcbioSingleCell)
library(splatter)
library(Seurat)
library(Matrix)
library(tidyverse)

# Check and make sure Python umap-learn is accessible to run UMAP.
# We're using this in the `Seurat::RunUMAP()` call below.
# Set `RETICULATE_PYTHON` to conda python binary in `~/.Renviron`.
# use_condaenv("steinbaugh")
stopifnot(py_module_available(module = "umap"))

organism <- "Homo sapiens"
release <- 92L

# splatter =====================================================================
# Use splatter to generate an example dataset with simulated counts.
# Note: These DE params are natural log scale.
params <- newSplatParams()
params <- setParam(params, "de.facLoc", 1)
params <- setParam(params, "de.facScale", .25)
params <- setParam(params, "dropout.type", "experiment")
params <- setParam(params, "dropout.mid", 3)
sce <- splatSimulate(params, group.prob = c(.5, .5), method = "groups")
# Prepare column data.
colData(sce) <- camel(colData(sce), rownames = TRUE, colnames = TRUE)
sce$cell <- NULL
# Add sampleID and sampleName columns.
sce$batch <- as.factor(camel(sce$batch))
sce$group <- as.factor(camel(sce$group))
sce$sampleID <- sce$group
sce$sampleName <- sce$sampleID
stopifnot("sampleName" %in% colnames(colData(sce)))
# Just keep raw counts.
assays(sce) <- assays(sce)["counts"]
# Prepare row data.
gr <- makeGRangesFromEnsembl(organism = organism, release = release)
rowRanges(sce) <- gr[seq_len(nrow(sce))]
# Sanitize to camel case.
colData(sce) <- camel(colData(sce))
metadata(sce) <- list()
# Set the interesting groups to sample name.
interestingGroups(sce) <- "sampleName"
sce <- metrics(sce, recalculate = TRUE)
sce <- filterCells(sce, minCellsPerGene = 25)

# seurat_small =================================================================
seurat_small <- sce %>%
    convertGenesToSymbols() %>%
    as("seurat") %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(resolution = seq(from = 0.4, to = 1.2, by = 0.4)) %>%
    RunTSNE(check_duplicates = FALSE) %>%
    RunUMAP() %>%
    SetAllIdent(id = "res.0.4")

# sce_small ====================================================================
# Convert rows (geneName) back to Ensembl IDs (geneID).
seurat_sce <- seurat_small %>%
    as("SingleCellExperiment") %>%
    convertSymbolsToGenes()
# Check that the conversion worked as expected.
stopifnot(identical(dimnames(sce), dimnames(seurat_sce)))
stopifnot(all(c("ident", "origIdent") %in% colnames(colData(seurat_sce))))
# Ensure that dimensional reduction data is slotted correctly.
stopifnot(identical(
    names(reducedDims(seurat_sce)),
    c("PCA", "TSNE", "UMAP")
))
# Overwrite the slots with the seurat data.
assays(sce) <- assays(seurat_sce)
colData(sce) <- colData(seurat_sce)
reducedDims(sce) <- reducedDims(seurat_sce)
# Remove genes with all zero counts prior to zinbwave.
sce <- filterCells(sce)
# Calculate zero weights using zinbwave (ZINB-WaVE negative binomial fit).
sce <- runZinbwave(sce)
sce_small <- sce

# all_markers_small ============================================================
all_markers_small <- SeuratMarkers(
    data = FindAllMarkers(seurat_small),
    GRanges = rowRanges(seurat_small)
)

# known_markers_small ==========================================================
g2s <- gene2symbol(seurat_small)
data <- tibble(
    cellType = as.factor(paste("cell_type", seq_len(2L), sep = "_")),
    geneID = head(g2s[["geneID"]], n = 2L),
    geneName = head(g2s[["geneName"]], n = 2L)
) %>%
    group_by(cellType)
known_markers_small <- new(
    Class = "CellTypeMarkers",
    as(data, "DataFrame"),
    metadata = list(
        version = packageVersion("pointillism"),
        organism = organism,
        ensemblRelease = release,
        date = Sys.Date()
    )
)

# Write out an example CSV that we can use to test `CellTypeMarkers()`.
export(
    x = known_markers_small,
    file = file.path("inst", "extdata", "cell_type_markers.csv")
)

# known_markers_detected_small =================================================
# FIXME Need to update this function
known_markers_detected_small <- knownMarkersDetected(
    all = all_markers_small,
    known = known_markers_small
)

# Save =========================================================================
devtools::use_data(
    sce_small,
    seurat_small,
    all_markers_small,
    known_markers_small,
    known_markers_detected_small,
    compress = "xz",
    overwrite = TRUE
)
