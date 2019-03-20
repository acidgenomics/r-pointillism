# Seurat example data
# 2019-03-20

library(pryr)
library(reticulate)
library(splatter)
library(Seurat)
library(Matrix)
library(tidyverse)
library(bcbioSingleCell)

# Check and make sure Python umap-learn is accessible to run UMAP.
# We're using this in the `Seurat::RunUMAP` call below.
# Set `RETICULATE_PYTHON` to conda python binary in `~/.Renviron`.
# This is not working consistently for me on Linux.

use_condaenv(
    condaenv = "reticulate",
    conda = "/usr/local/miniconda3/bin/conda"
)

# umap-learn via reticulate doesn't work well with conda.
# https://github.com/satijalab/seurat/issues/486

stopifnot(
    identical(basename(Sys.getenv("RETICULATE_PYTHON")), "python"),
    py_module_available(module = "umap")
)

# # Restrict object size to 1 MB.
# Use `pryr::object_size` instead of `utils::object.size`.
limit <- structure(1e6, class = "object_size")

data(pbmc_small, package = "Seurat")
object_size(pbmc_small)
stopifnot(object_size(pbmc_small) < limit)

# seurat_small =================================================================
seurat_small <- pbmc_small
# Note that this step requires umap-learn via reticulate.
seurat_small <- RunUMAP(seurat_small)
object_size(seurat_small)
stopifnot(object_size(seurat_small) < limit)
validObject(seurat_small)

# `Seurat::pbmc_small` gene symbols map to GRCh37.
gr <- makeGRangesFromEnsembl("Homo sapiens", genomeBuild = "GRCh37")
x <- rownames(seurat_small)
table <- gr$geneName %>%
    as.character() %>%
    make.unique()
names(gr) <- table
stopifnot(all(x %in% table))
which <- match(x = x, table = table)
gr <- gr[which]
rowRanges(seurat_small) <- gr

# all_markers_small ============================================================
markers <- FindAllMarkers(seurat_small)
ranges <- rowRanges(seurat_small)
all_markers_small <- SeuratMarkersPerCluster(object = markers, ranges = ranges)

# known_markers_small ==========================================================
known_markers_small <- KnownMarkers(
    markers = all_markers_small,
    known = cell_type_markers$homoSapiens
)
export(
    x = known_markers_small,
    file = file.path("inst", "extdata", "cell_type_markers.csv")
)

# Save =========================================================================
usethis::use_data(
    seurat_small,
    all_markers_small,
    known_markers_small,
    compress = "xz",
    overwrite = TRUE
)
