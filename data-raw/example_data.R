# SingleCellExperiment Example Data
# 2018-10-17

# Check and make sure Python umap-learn is accessible to run UMAP.
# We're using this in the `Seurat::RunUMAP()` call below.
# Set `RETICULATE_PYTHON` to conda python binary in `~/.Renviron`.
library(reticulate)
# use_condaenv("bioinfo")
# Unable to find conda binary. Is Anaconda installed?
stopifnot(py_module_available(module = "umap"))

library(splatter)
library(Seurat)
library(Matrix)
library(tidyverse)

library(bcbioSingleCell)
sce <- bcbioSingleCell::cellranger_small
organism <- metadata(sce)$organism
release <- 84

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

# all_markers_small ============================================================
all_markers_small <- SeuratMarkers(
    markers = FindAllMarkers(seurat_small),
    ranges = rowRanges(seurat_small)
)

# known_markers_small ==========================================================
all <- all_markers_small
# FIXME Need to make `CellTypeMarkers()` a generic.
known <- new(
    Class = "CellTypeMarkers",
    DataFrame(
        cellType = as.factor(paste("cell_type", seq_len(2L), sep = "_")),
        geneID = as.character(head(all$ranges$geneID, n = 2L)),
        geneName = as.character(head(all$ranges$geneName, n = 2L))
    ),
    metadata = list(
        version = packageVersion("pointillism"),
        organism = organism,
        ensemblRelease = release,
        date = Sys.Date()
    )
)
# Write out an example CSV that we can use to test `CellTypeMarkers()`.
export(
    x = known,
    file = file.path("inst", "extdata", "cell_type_markers.csv")
)

known_markers_small <- knownMarkers(all = all, known = known)

# Save =========================================================================
usethis::use_data(
    seurat_small,
    all_markers_small,
    known_markers_small,
    compress = "xz",
    overwrite = TRUE
)
