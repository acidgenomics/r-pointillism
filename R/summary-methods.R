#' @name summary
#' @importFrom S4Vectors summary
#' @inherit basejump::summary
#' @note Updated 2020-01-30.
#'
#' @inheritParams acidroxygen::params
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



#' @rdname summary
#' @name summary
#' @importFrom S4Vectors summary
#' @usage summary(object, ...)
#' @export
NULL



## Updated 2019-07-31.
`summary,SeuratMarkers` <-  # nolint
    function(object) {
        ## Metadata.
        m <- metadata(object)
        ## Row ranges metadata.
        rrm <- metadata(object[["ranges"]])
        showSlotInfo(list(
            alpha = m[["alpha"]],
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



## Updated 2019-07-31.
`summary,SeuratMarkersPerCluster` <-  # nolint
    function(object) {
        cat(paste(length(object), "clusters"), sep = "\n")
        data <- do.call(what = rbind, args = object)
        ## Metadata.
        m <- metadata(data)
        ## Row ranges metadata.
        rrm <- metadata(data[["ranges"]])
        showSlotInfo(list(
            alpha = m[["alpha"]],
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
