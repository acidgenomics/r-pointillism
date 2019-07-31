## Seurat example markers return data.
## Updated 2019-07-31.

library(usethis)
library(Seurat)  # 3.0.2

data(cellTypeMarkersList)
data(Seurat, package = "acidtest")

## seurat_all_markers ===========================================================
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

## known_markers_small ==========================================================
seuratKnownMarkers <- KnownMarkers(
    markers = seuratAllMarkers,
    known = cellTypeMarkersList[["homoSapiens"]]
)

## Save =========================================================================
use_data(
    seuratAllMarkers,
    seuratKnownMarkers,
    compress = "xz",
    overwrite = TRUE
)
