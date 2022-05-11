tmpenv <- new.env()
data(seurat, smpc, envir = tmpenv)
data(KnownMarkers, package = "AcidTest", envir = tmpenv)
objs <- list()
objs[["KnownMarkers"]] <- get("KnownMarkers", envir = tmpenv)
objs[["Seurat"]] <- get("seurat", envir = tmpenv)
objs[["SeuratMarkersPerCluster"]] <- get("smpc", envir = tmpenv)
objs[["SingleCellExperiment"]] <- as(objs[["Seurat"]], "SingleCellExperiment")
objs[["genes"]] <- head(rownames(objs[["Seurat"]]))
rm(tmpenv)

## nolint start
CellTypeMarkers <- AcidSingleCell::CellTypeMarkers
## nolint end
