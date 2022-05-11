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

SingleCellExperiment <- as(Seurat, "SingleCellExperiment")
genes <- head(rownames(Seurat))

## FIXME Rework this...take out.
objects <- list(
    "SingleCellExperiment" = SingleCellExperiment,
    "Seurat" = Seurat
)

## nolint start
CellTypeMarkers <- AcidSingleCell::CellTypeMarkers
## nolint end
