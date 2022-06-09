tmpenv <- new.env()
data(seurat, smpc, envir = tmpenv)
data(
    KnownMarkers,
    Seurat,
    SeuratMarkersPerCluster,
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
objs[["genes"]] <-
    head(rownames(objs[["Seurat"]]))
rm(tmpenv)

## nolint start
CellTypeMarkers <- AcidSingleCell::CellTypeMarkers
## nolint end
