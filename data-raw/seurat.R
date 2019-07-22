## Seurat example data
## Updated 2019-07-16.

library(pryr)
library(reticulate)
library(Matrix)
library(splatter)
library(Seurat)  # 3.0
library(bcbioSingleCell)
library(tidyverse)

data(seurat, package = "acidtest")

## seurat_all_markers ===========================================================
seurat_all_markers <- SeuratMarkersPerCluster(
    object = FindAllMarkers(seurat),
    ranges = rowRanges(seurat)
)

## known_markers_small ==========================================================
data(cell_type_markers)
seurat_known_markers <- KnownMarkers(
    markers = seurat_all_markers,
    known = cell_type_markers$homoSapiens
)
export(
    x = seurat_known_markers,
    file = file.path("inst", "extdata", "cell_type_markers.csv")
)

## Save =========================================================================
usethis::use_data(
    seurat,
    seurat_all_markers,
    seurat_known_markers,
    compress = "xz",
    overwrite = TRUE
)
