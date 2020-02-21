## nocov start



#' @name deprecated
#' @inherit basejump::deprecated
#' @keywords internal
NULL



## v0.4.0 ======================================================================
#' @rdname deprecated
#' @export
clusterID <- function(...) {
    .Deprecated("clusters")
    clusters(...)
}



## v0.4.5 ======================================================================
#' @importFrom acidgenerics plotDot
NULL

`plotDot,ANY` <-  # nolint
    function(object, ...) {
        .Deprecated("plotDots")
        plotDots(object, ...)
    }

#' @rdname deprecated
#' @export
setMethod(
    f = "plotDot",
    signature = signature("ANY"),
    definition = `plotDot,ANY`
)



## nocov end
