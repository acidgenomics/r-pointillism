#' @name topMarkers
#' @inherit AcidGenerics::topMarkers
#' @note Updated 2020-01-30.
#'
#' @inheritParams AcidRoxygen::params
#' @param n `integer(1)`.
#'   Number of genes per cluster.
#' @param direction `character(1)`.
#'   Whether to include upregulated (`"up"`; positive LFC), downregulated
#'   (`"down"`; negative LFC) or `"both"` directions of association per cluster.
#' @param ... Additional arguments.
#'
#' @examples
#' data(seurat_all_markers)
#'
#' ## SeuratMarkersPerCluster ====
#' object <- seurat_all_markers
#' x <- topMarkers(object, n = 2L)
#' print(x)
NULL



## Note that the validity method checks for sorting by adjusted P value.
## Updated 2020-02-21.
`topMarkers,SeuratMarkersPerCluster` <-  # nolint
    function(
        object,
        n = 10L,
        direction
    ) {
        validObject(object)
        assert(isInt(n))
        direction <- match.arg(direction)
        cli_alert_info(switch(
            EXPR = direction,
            "both" = "Including both up- and down-regulated markers.",
            "down" = "Including downregulated markers.",
            "up" = "Including upregulated markers."
        ))
        x <- object
        x <- SplitDataFrameList(lapply(
            X = x,
            FUN = function(x) {
                ## Subset to positive or negative correlation, if desired.
                if (identical(direction, "up")) {
                    keep <- x[["avgLogFC"]] > 0L
                    x <- x[keep, , drop = FALSE]
                } else if (identical(direction, "down")) {
                    keep <- x[["avgLogFC"]] < 0L
                    x <- x[keep, , drop = FALSE]
                }
                x <- head(x, n = n)
                x
            }
        ))
        x <- unlist(x, recursive = FALSE, use.names = FALSE)
        ranges <- x[["ranges"]]
        x[["ranges"]] <- NULL
        x[["geneID"]] <- as.character(mcols(ranges)[["geneID"]])
        x[["geneName"]] <- as.character(mcols(ranges)[["geneName"]])
        x <- mutateIf(x, is.character, as.factor)
        x
    }

args <- "direction"
formals(`topMarkers,SeuratMarkersPerCluster`)[args] <-
    .formalsList[args]
rm(args)



#' @rdname topMarkers
#' @export
setMethod(
    f = "topMarkers",
    signature = signature("SeuratMarkersPerCluster"),
    definition = `topMarkers,SeuratMarkersPerCluster`
)
