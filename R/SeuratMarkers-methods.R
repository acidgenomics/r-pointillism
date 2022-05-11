#' Sanitize Seurat markers
#'
#' This generator function is designed to take the original return from a Seurat
#' marker analysis and add corresponding gene annotations.
#'
#' @name SeuratMarkers
#' @note Updated 2022-05-11.
#'
#' @details
#' For [Seurat::FindAllMarkers()] return, rownames are correctly returned
#' in the `gene` column.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @param object
#' Unmodified Seurat marker return `data.frame`.
#' - `SeuratMarkers()`: [Seurat::FindMarkers()].
#' - `SeuratMarkersPerCluster()`: [Seurat::FindAllMarkers()].
#'
#' @param ranges `GenomicRanges`.
#' Gene annotations. Names must correspond to the rownames. The function will
#' automatically subset the ranges and arrange them alphabetically.
#'
#' @return `SeuratMarkers`.
#'
#' @examples
#' data(seurat)
#'
#' ## Seurat ====
#' object <- seurat
#' ranges <- rowRanges(object)
#'
#' ## `FindMarkers()` return.
#' invisible(capture.output({
#'     markers <- Seurat::FindMarkers(
#'         object = object,
#'         ident.1 = "1",
#'         ident.2 = NULL
#'     )
#' }))
#' x <- SeuratMarkers(object = markers, ranges = ranges)
#' summary(x)
#'
#' ## `FindAllMarkers()` return.
#' invisible(capture.output(suppressWarnings({
#'     markers <- Seurat::FindAllMarkers(object)
#' })))
#' x <- SeuratMarkersPerCluster(object = markers, ranges = ranges)
#' summary(x)
NULL



## Updated 2021-03-03.
`SeuratMarkers,data.frame` <- # nolint
    function(object,
             ranges,
             alphaThreshold = 0.05) {
        assert(
            hasRows(object),
            hasRownames(object),
            is(ranges, "GenomicRanges"),
            isSubset(
                x = c("geneId", "geneName"),
                y = colnames(mcols(ranges))
            ),
            isAlpha(alphaThreshold)
        )
        x <- as(object, "DataFrame")
        ## Update legacy columns.
        if (isSubset("avg_diff", colnames(x))) {
            alertWarning(sprintf(
                paste(
                    "Renaming legacy {.var %s} column to {.var %s}",
                    "(changed in {.pkg %s} v%s)."
                ),
                "avg_diff", "avg_log2FC", "Seurat", "2.1"
            ))
            colnames(x)[colnames(x) == "avg_diff"] <- "avg_log2FC"
        }
        if (isSubset("avg_logFC", colnames(x))) {
            alertWarning(sprintf(
                paste(
                    "Renaming legacy {.var %s} column to {.var %s}",
                    "(changed in {.pkg %s} v%s)."
                ),
                "avg_logFC", "avg_log2FC", "Seurat", "4.0"
            ))
            colnames(x)[colnames(x) == "avg_logFC"] <- "avg_log2FC"
        }
        ## Sanitize the column names in to camel case.
        cols <- snakeCase(colnames(x))
        ## Handle "avg_log2FC" column.
        cols <- gsub(
            pattern = "([0-9]+)([a-z]+)",
            replacement = "\\1_\\2",
            x = cols
        )
        cols <- camelCase(cols, strict = TRUE)
        colnames(x) <- cols
        fun <- ifelse(
            test = isSubset(c("cluster", "gene"), colnames(x)),
            yes = "FindAllMarkers",
            no = "FindMarkers"
        )
        alertInfo(sprintf("{.fun %s} return detected.", fun))
        ## Map the Seurat matrix rownames to `rownames` column in tibble.
        ## NOTE Using "name" here instead of "geneName", which is intended
        ## to map to the rownames, which can be altered by `make.names`.
        switch(
            EXPR = fun,
            "FindMarkers" = {
                perCluster <- FALSE
                x[["name"]] <- rownames(x)
            },
            "FindAllMarkers" = {
                perCluster <- TRUE
                colnames(x)[colnames(x) == "gene"] <- "name"
            }
        )
        rownames(x) <- NULL
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
            "avgLog2Fc",
            "padj",
            "pvalue"
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
        metadata(x) <- append(
            x = .prototypeMetadata,
            values = list(
                "alphaThreshold" = alphaThreshold,
                "sessionInfo" = sessionInfo()
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



## Updated 2019-08-06.
`SeuratMarkersPerCluster,data.frame` <- # nolint
    `SeuratMarkers,data.frame`



#' @rdname SeuratMarkers
#' @export
setMethod(
    f = "SeuratMarkers",
    signature = signature(object = "data.frame"),
    definition = `SeuratMarkers,data.frame`
)



#' @rdname SeuratMarkers
#' @export
setMethod(
    f = "SeuratMarkersPerCluster",
    signature = signature(object = "data.frame"),
    definition = `SeuratMarkersPerCluster,data.frame`
)
