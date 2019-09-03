## FIXME Update bioverbs return.
## FIXME Need to generate a multi-sample example.



#' @name clusterCellCountsPerSample
#' @inherit bioverbs::clusterCellCountsPerSample
#' @note Updated 2019-09-03.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(
#'     Seurat,
#'     cell_data_set,
#'     package = "acidtest"
#' )
#'
#' ## Seurat ====
#' object <- Seurat
#' colData(object)[["sampleID"]] <- as.factor(paste0("sample", seq(2L)))
#' colData(object)[["sampleName"]] <- colData(object)[["sampleID"]]
#' x <- clusterCellCountsPerSample(object)
#' print(x)
#'
#' ## cell_data_set ====
#' object <- cell_data_set
#' x <- clusterCellCountsPerSample(object)
#' print(x)
NULL



#' @rdname clusterCellCountsPerSample
#' @name clusterCellCountsPerSample
#' @importFrom bioverbs clusterCellCountsPerSample
#' @usage clusterCellCountsPerSample(object, ...)
#' @export
NULL



## nolint start
## dplyr approach:
## > x <- metrics(object)
## > x <- arrange(x, !!!syms(cols))
## > x <- group_by(x, !!!syms(cols))
## > x <- summarize(x, n = n())
## > x <- ungroup(x)
## > x <- arrange(x, !!!syms(cols))
## > x <- group_by(x, !!sym("sampleName"))
## > x <- mutate(x, ratio = !!sym("n") / sum(!!sym("n")))
## nolint end



## Updated 2019-08-03.
`clusterCellCountsPerSample,SingleCellExperiment` <-  # nolint
    function(object) {
        assert(.hasClusters(object))
        x <- metrics(object, return = "DataFrame")
        cols <- c("ident", "sampleName")
        assert(isSubset(cols, colnames(x)))
        x <- x[, cols, drop = FALSE]
        tbl <- table(x)
        ## Get the number of cells per ident.
        ## We'll use this in a join below to calculate ratio.
        nPerIdent <- rowSums(tbl)
        nPerIdent <- DataFrame(
            ident = names(nPerIdent),
            nPerIdent = nPerIdent
        )
        x <- melt(tbl, colnames = c(cols, "n"))
        x <- leftJoin(x, nPerIdent, by = "ident")
        x[["ratio"]] <- x[["n"]] / x[["nPerIdent"]]
        x <- x[order(x[["ident"]], x[["sampleName"]]), , drop = FALSE]
        x
    }



#' @rdname clusterCellCountsPerSample
#' @export
setMethod(
    f = "clusterCellCountsPerSample",
    signature = signature("SingleCellExperiment"),
    definition = `clusterCellCountsPerSample,SingleCellExperiment`
)



## Updated 2019-07-31.
`clusterCellCountsPerSample,Seurat` <-  # nolint
    `clusterCellCountsPerSample,SingleCellExperiment`



#' @rdname clusterCellCountsPerSample
#' @export
setMethod(
    f = "clusterCellCountsPerSample",
    signature = signature("Seurat"),
    definition = `clusterCellCountsPerSample,Seurat`
)



## Updated 2019-08-05.
`clusterCellCountsPerSample,cell_data_set` <-  # nolint
    `clusterCellCountsPerSample,SingleCellExperiment`



#' @rdname clusterCellCountsPerSample
#' @export
setMethod(
    f = "clusterCellCountsPerSample",
    signature = signature("cell_data_set"),
    definition = `clusterCellCountsPerSample,cell_data_set`
)
