## Updated 2019-07-31.
.prototypeMetadata <- list(
    version = packageVersion(packageName()),
    date = Sys.Date()
)



#' Cell-cycle markers
#'
#' Data provenence information, including the organism and Ensembl release are
#' defined in [`metadata()`][S4Vectors::metadata].
#'
#' @export
#' @note Updated 2019-09-03.
#'
#' @return `CellTypeMarkers`
setClass(
    Class = "CellCycleMarkers",
    contains = "CompressedSplitDataFrameList"
)
setValidity(
    Class = "CellCycleMarkers",
    method = function(object) {
        validate(
            areSetEqual(
                x = colnames(object[[1L]]),
                y = c("geneId", "geneName", "phase")
            ),
            isSubset(
                x = c("version", "organism", "ensemblRelease", "date"),
                y = names(metadata(object))
            )
        )
    }
)



#' Cell-type markers
#'
#' Data provenence information, including the organism and Ensembl release are
#' defined in [`metadata()`][S4Vectors::metadata].
#'
#' @export
#' @note Updated 2019-09-03.
#'
#' @return `CellTypeMarkers`
setClass(
    Class = "CellTypeMarkers",
    contains = "CompressedSplitDataFrameList"
)
setValidity(
    Class = "CellTypeMarkers",
    method = function(object) {
        validate(
            areSetEqual(
                x = colnames(object[[1L]]),
                y = c("cellType", "geneId", "geneName")
            ),
            isSubset(
                x = c("date", "ensemblRelease", "organism", "version"),
                y = names(metadata(object))
            )
        )
    }
)



#' Known markers
#'
#' Class containing known markers detected.
#'
#' Results are grouped by `cellType` column and arranged by adjusted *P* value
#' (`padj`).
#'
#' @export
#' @note Updated 2020-10-12.
#'
#' @return `KnownMarkers`.
setClass(
    Class = "KnownMarkers",
    contains = "DataFrame"
)
## Consider requiring "avgLogFC" column.
setValidity(
    Class = "KnownMarkers",
    method = function(object) {
        validate(
            isSubset(
                x = c(
                    "cellType",
                    "cluster",
                    "geneId",
                    "geneName",
                    "name",
                    "padj",
                    "pvalue"
                ),
                y = colnames(object)
            ),
            isSubset(
                x = c("alphaThreshold", "date", "version"),
                y = names(metadata(object))
            )
        )
    }
)



#' Seurat markers
#'
#' Class containing essential elements for marker gene analysis.
#'
#' Results are arranged by adjusted *P* value (`padj`).
#'
#' @export
#' @note Updated 2019-09-03.
#'
#' @return `SeuratMarkers`.
setClass(
    Class = "SeuratMarkers",
    contains = "DataFrame"
)
setValidity(
    Class = "SeuratMarkers",
    method = function(object) {
        validate(
            hasRownames(object),
            identical(
                x = sort(colnames(object)),
                y = c(
                    "avgLogFC",
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
#' @note Updated 2019-09-03.
#'
#' @return `SeuratMarkersPerCluster`.
setClass(
    Class = "SeuratMarkersPerCluster",
    contains = "CompressedSplitDataFrameList"
)
setValidity(
    Class = "SeuratMarkersPerCluster",
    method = function(object) {
        validate(
            all(grepl("^cluster", names(object))),
            identical(
                x = sort(colnames(object)[[1L]]),
                y = c(
                    "avgLogFC",
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
