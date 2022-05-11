## Fix for edgeR partial match warning.
options(
    "warnPartialMatchAttr" = FALSE,
    "warnPartialMatchDollar" = FALSE
)

tmpenv <- new.env()
data(
    ## > cell_data_set,
    Seurat,
    SeuratMarkersPerCluster,
    envir = tmpenv
)
data(
    KnownMarkers,
    package = "AcidTest",
    envir = tmpenv
)
objs <- list()
objs[["KnownMarkers"]] <-
    get("KnownMarkers", envir = tmpenv)
objs[["Seurat"]] <-
    get("Seurat", envir = tmpenv)
objs[["SeuratMarkersPerCluster"]] <-
    get("SeuratMarkersPerCluster", envir = tmpenv)
objs[["SingleCellExperiment"]] <-
    as(objs[["Seurat"]], "SingleCellExperiment")
objs[["genes"]] <- head(rownames(objs[["Seurat"]]))
rm(tmpenv)

## nolint start
CellTypeMarkers <- AcidSingleCell::CellTypeMarkers
## nolint end
