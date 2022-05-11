## Fix for edgeR partial match warning.
options(
    "warnPartialMatchAttr" = FALSE,
    "warnPartialMatchDollar" = FALSE
)

data(
    ## > cell_data_set,
    Seurat,
    SeuratMarkersPerCluster,
    envir = environment()
)
data(
    KnownMarkers,
    package = "AcidTest",
    envir = environment()
)

sce <- as(Seurat, "SingleCellExperiment")
genes <- head(rownames(Seurat))
seurat <- Seurat
rm(Seurat)
smpc <- SeuratMarkersPerCluster
rm(SeuratMarkersPerCluster)

## nolint start
CellTypeMarkers <- AcidSingleCell::CellTypeMarkers
## nolint end
