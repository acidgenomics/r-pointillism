#' @name summary
#' @importFrom S4Vectors summary
#' @inherit basejump::summary
#' @inheritParams basejump::params
#' @examples
#' data(rse, txse, package = "acidtest")
#'
#' ## Gene2Symbol ====
#' x <- Gene2Symbol(rse)
#' summary(x)
#'
#' ## Tx2Gene ====
#' x <- Tx2Gene(txse)
#' summary(x)
NULL



#' @rdname summary
#' @name summary
#' @importFrom S4Vectors summary
#' @export
NULL



summary.SeuratMarkers <-  # nolint
    function(object) {
        # Metadata.
        m <- metadata(object)
        # Row ranges metadata.
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
    definition = summary.SeuratMarkers
)



summary.SeuratMarkersPerCluster <-  # nolint
    function(object) {
        cat(paste(length(object), "clusters"), sep = "\n")
        data <- do.call(what = rbind, args = object)
        # Metadata.
        m <- metadata(data)
        # Row ranges metadata.
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
    definition = summary.SeuratMarkersPerCluster
)
