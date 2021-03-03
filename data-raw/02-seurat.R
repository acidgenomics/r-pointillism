## Seurat example markers return data.
## Updated 2021-03-03.

library(usethis)
library(Seurat)  # 4.0.0

data(cellTypeMarkersList)
data(Seurat, package = "AcidTest")

seuratAllMarkers <- withCallingHandlers(
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

seuratKnownMarkers <- KnownMarkers(
    markers = seuratAllMarkers,
    known = cellTypeMarkersList[["homoSapiens"]]
)

use_data(
    seuratAllMarkers,
    seuratKnownMarkers,
    compress = "xz",
    overwrite = TRUE
)
