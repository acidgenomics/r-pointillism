#' @name summary
#' @inherit base::summary description return title
#' @note Updated 2020-10-12.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(seurat_all_markers)
#'
#' ## SeuratMarkersPerCluster ====
#' object <- seurat_all_markers
#' class(object)
#' summary(object)
NULL



## Updated 2020-10-12.
`summary,SeuratMarkers` <-  # nolint
    function(object) {
        ## Metadata.
        m <- metadata(object)
        ## Row ranges metadata.
        rrm <- metadata(object[["ranges"]])
        showSlotInfo(list(
            alphaThreshold = m[["alphaThreshold"]],
            organism = rrm[["organism"]],
            genomeBuild = rrm[["genomeBuild"]],
            ensemblRelease = rrm[["ensemblRelease"]],
            version = as.character(m[["version"]]),
            date = m[["date"]]
        ))
    }



#' @rdname summary
#' @export
setMethod(
    f = "summary",
    signature = signature("SeuratMarkers"),
    definition = `summary,SeuratMarkers`
)



## Updated 2020-10-12.
`summary,SeuratMarkersPerCluster` <-  # nolint
    function(object) {
        cat(paste(length(object), "clusters"), sep = "\n")
        data <- do.call(what = rbind, args = object)
        ## Metadata.
        m <- metadata(data)
        ## Row ranges metadata.
        rrm <- metadata(data[["ranges"]])
        showSlotInfo(list(
            alphaThreshold = m[["alphaThreshold"]],
            organism = rrm[["organism"]],
            genomeBuild = rrm[["genomeBuild"]],
            ensemblRelease = rrm[["ensemblRelease"]],
            version = as.character(m[["version"]]),
            date = m[["date"]]
        ))
    }



#' @rdname summary
#' @export
setMethod(
    f = "summary",
    signature = signature("SeuratMarkersPerCluster"),
    definition = `summary,SeuratMarkersPerCluster`
)
