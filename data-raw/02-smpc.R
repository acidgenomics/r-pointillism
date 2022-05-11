suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(AcidSingleCell) # 0.3.0
    library(Seurat) # 4.1.0
})
load_all() # helpers = FALSE
data(seurat)
data(cellTypeMarkersList, package = "AcidSingleCell")
object <- seurat
smpc <-
    withCallingHandlers(
        expr = SeuratMarkersPerCluster(
            object = FindAllMarkers(object),
            ranges = rowRanges(object)
        ),
        warning = function(w) {
            if (grepl("cannot compute exact p-value with ties", w)) {
                invokeRestart("muffleWarning")
            } else {
                w
            }
        }
    )
use_data(smpc, compress = "xz", overwrite = TRUE)
