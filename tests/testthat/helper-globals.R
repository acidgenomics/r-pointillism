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

objs <- list()
objs[["KnownMarkers"]] <- KnownMarkers
objs[["Seurat"]] <- Seurat
objs[["SeuratMarkersPerCluster"]] <- SeuratMarkersPerCluster
objs[["SingleCellExperiment"]] <-
    as(objs[["Seurat"]], "SingleCellExperiment")
objs[["genes"]] <- head(rownames(objs[["Seurat"]]))

rm(
    KnownMarkers,
    Seurat,
    SeuratMarkersPerCluster
)

## nolint start
CellTypeMarkers <- AcidSingleCell::CellTypeMarkers
## nolint end
