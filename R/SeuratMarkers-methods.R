#' Sanitize Seurat markers
#'
#' This generator function is designed to take the original return from a Seurat
#' marker analysis and add corresponding gene annotations.
#'
#' @name SeuratMarkers
#' @note For [Seurat::FindAllMarkers()] return, rownames are correctly returned
#'   in the `gene` column.
#' @note Updated 2020-01-30.
#'
#' @inheritParams acidroxygen::params
#'
#' @param object
#'   Unmodified Seurat marker return `data.frame`.
#'   - `SeuratMarkers()`: [Seurat::FindMarkers()].
#'   - `SeuratMarkersPerCluster()`: [Seurat::FindAllMarkers()].
#' @param ranges `GRanges`.
#'   Gene annotations. Names must correspond to the rownames. The function will
#'   automatically subset the ranges and arrange them alphabetically.
#'
#' @return `SeuratMarkers`.
#'
#' @examples
#' data(Seurat, package = "acidtest")
#'
#' ## Seurat ====
#' object <- Seurat
#' ranges <- rowRanges(object)
#'
#' ## `FindMarkers()` return.
#' invisible(capture.output(
#'     markers <- Seurat::FindMarkers(
#'         object = object,
#'         ident.1 = "1",
#'         ident.2 = NULL
#'     )
#' ))
#' x <- SeuratMarkers(object = markers, ranges = ranges)
#' summary(x)
#'
#' ## `FindAllMarkers()` return.
#' invisible(capture.output(suppressWarnings(
#'     markers <- Seurat::FindAllMarkers(object)
#' )))
#' x <- SeuratMarkersPerCluster(object = markers, ranges = ranges)
#' summary(x)
NULL



## Updated 2020-01-30.
`SeuratMarkers,data.frame` <-  # nolint
    function(
        object,
        ranges,
        alpha = 0.05
    ) {
        assert(
            hasRows(object),
            hasRownames(object),
            is(ranges, "GRanges"),
            isSubset(
                x = c("geneID", "geneName"),
                y = colnames(mcols(ranges))
            ),
            isAlpha(alpha)
        )
        ## Detect function from column names.
        cols <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")
        if (identical(colnames(object), cols)) {
            perCluster <- FALSE
            fun <- "FindMarkers"
        } else if (identical(colnames(object), c(cols, "cluster", "gene"))) {
            perCluster <- TRUE
            fun <- "FindAllMarkers"
        }
        cli_alert_info(sprintf("{.fun %s} return detected.", fun))
        ## Sanitize markers.
        x <- as(object, "DataFrame")
        x <- camelCase(x)
        ## Map the Seurat matrix rownames to `rownames` column in tibble.
        if (identical(fun, "FindMarkers")) {
            x[["name"]] <- rownames(x)
        } else if (identical(fun, "FindAllMarkers")) {
            colnames(x)[colnames(x) == "gene"] <- "name"
        }
        rownames(x) <- NULL
        ## Update legacy columns.
        if (isSubset("avgDiff", colnames(x))) {
            cli_alert_warning(paste(
                "Renaming legacy {.var avgDiff} column to {.var avgLogFC}",
                "(changed in {.pkg Seurat} v2.1)."
            ))
            colnames(x)[colnames(x) == "avgDiff"] <- "avgLogFC"
        }
        ## Rename P value columns to match DESeq2 conventions.
        if (isSubset("pVal", colnames(x))) {
            colnames(x)[colnames(x) == "pVal"] <- "pvalue"
        }
        if (isSubset("pValAdj", colnames(x))) {
            colnames(x)[colnames(x) == "pValAdj"] <- "padj"
        }
        ## Ensure that required columns are present.
        requiredCols <- c(
            "name",
            "pct1",
            "pct2",
            "avgLogFC",  # Seurat v2.1.
            "padj",
            "pvalue"     # Renamed from `p_val`.
        )
        assert(isSubset(requiredCols, colnames(x)))
        ## Bind ranges as column.
        assert(isSubset(unique(x[["name"]]), names(ranges)))
        x[["ranges"]] <- ranges[x[["name"]]]
        ## Arrange by adjusted P value.
        x <- x[order(x[["padj"]]), sort(colnames(x)), drop = FALSE]
        ## Split by cluster, if applicable.
        if (isTRUE(perCluster)) {
            x <- split(x, f = x[["cluster"]])
        }
        ## Add metadata and return.
        metadata(x) <- c(
            .prototypeMetadata,
            list(
                alpha = alpha,
                sessionInfo = session_info(include_base = TRUE)
            )
        )
        if (isTRUE(perCluster)) {
            names(x) <- paste0("cluster", names(x))
            new(Class = "SeuratMarkersPerCluster", x)
        } else {
            rownames(x) <- x[["name"]]
            x[["name"]] <- NULL
            new(Class = "SeuratMarkers", x)
        }
    }



#' @rdname SeuratMarkers
#' @export
setMethod(
    f = "SeuratMarkers",
    signature = signature("data.frame"),
    definition = `SeuratMarkers,data.frame`
)



## Updated 2019-08-06.
`SeuratMarkersPerCluster,data.frame` <-  # nolint
    `SeuratMarkers,data.frame`



#' @rdname SeuratMarkers
#' @export
setMethod(
    f = "SeuratMarkersPerCluster",
    signature = signature("data.frame"),
    definition = `SeuratMarkersPerCluster,data.frame`
)
