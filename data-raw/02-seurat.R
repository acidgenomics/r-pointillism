## Seurat example markers return data.
## Updated 2020-10-12.

library(usethis)
library(Seurat)   # 4.0.0

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
    markers = seurat_all_markers,
    known = cell_type_markers_list[["homoSapiens"]]
)

use_data(
    seurat_all_markers,
    seurat_known_markers,
    compress = "xz",
    overwrite = TRUE
)
