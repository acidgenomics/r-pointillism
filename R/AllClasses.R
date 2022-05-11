#' Seurat markers
#'
#' Class containing essential elements for marker gene analysis.
#'
#' Results are arranged by adjusted *P* value (`padj`).
#'
#' @export
#' @note Updated 2022-05-11.
#'
#' @return `SeuratMarkers`.
setClass(
    Class = "SeuratMarkers",
    contains = "DFrame"
)
setValidity(
    Class = "SeuratMarkers",
    method = function(object) {
        validate(
            hasRownames(object),
            identical(
                x = sort(colnames(object)),
                y = c(
                    "avgLog2Fc",
                    "padj",
                    "pct1",
                    "pct2",
                    "pvalue",
                    "ranges"
                )
            ),
            ## Ensure sorting by adjusted P value.
            identical(object[["padj"]], sort(object[["padj"]]))
        )
    }
)



#' Seurat markers per cluster
#'
#' Class containing essential elements for marker per cluster analysis.
#'
#' Results are split per `cluster` and arranged by adjusted *P* value (`padj`).
#'
#' @export
#' @note Updated 2022-05-11.
#'
#' @return `SeuratMarkersPerCluster`.
setClass(
    Class = "SeuratMarkersPerCluster",
    contains = "CompressedSplitDFrameList"
)
setValidity(
    Class = "SeuratMarkersPerCluster",
    method = function(object) {
        validate(
            all(grepl("^cluster", names(object))),
            identical(
                x = sort(colnames(object)[[1L]]),
                y = c(
                    "avgLog2Fc",
                    "cluster",
                    "name",
                    "padj",
                    "pct1",
                    "pct2",
                    "pvalue",
                    "ranges"
                )
            ),
            ## Ensure sorting by adjusted P value.
            all(bapply(
                X = object,
                FUN = function(x) {
                    identical(x[["padj"]], sort(x[["padj"]]))
                }
            ))
        )
    }
)
