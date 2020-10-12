## Seurat example markers return data.
## Updated 2020-07-24.

library(usethis)
library(Seurat)   # 3.2.0

data(cell_type_markers_list)
data(Seurat, package = "AcidTest")

## seurat_all_markers ===========================================================
seurat_all_markers <- withCallingHandlers(
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
seurat_known_markers <- KnownMarkers(
    markers = seurat_all_markers,
    known = cell_type_markers_list[["homoSapiens"]]
)

## Save =========================================================================
use_data(
    seurat_all_markers,
    seurat_known_markers,
    compress = "xz",
    overwrite = TRUE
)
