suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(AcidSingleCell) # 0.3.0
    library(Seurat) # 4.1.0
})
load_all()
data(cellTypeMarkersList, package = "AcidSingleCell")
data(Seurat)
SeuratMarkersPerCluster <-  # nolint
    withCallingHandlers(
        expr = SeuratMarkersPerCluster(
            object = FindAllMarkers(Seurat),
            ranges = rowRanges(Seurat)
        ),
        warning = function(w) {
            if (grepl("cannot compute exact p-value with ties", w)) {
                invokeRestart("muffleWarning")
            } else {
                w
            }
        }
    )
use_data(
    SeuratMarkersPerCluster,
    compress = "xz",
    overwrite = TRUE
)
