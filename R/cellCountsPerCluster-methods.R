## FIXME Need to improve working example to include 2 samples.



#' @name cellCountsPerCluster
#' @inherit bioverbs::cellCountsPerCluster
#' @note Updated 2019-09-01.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @return `DataFrame`.
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
#' ## Simulate a multi-sample dataset.
#' colData(object)[["sampleID"]] <- as.factor(paste0("sample", seq(2L)))
#' colData(object)[["sampleName"]] <- colData(object)[["sampleID"]]
#' x <- cellCountsPerCluster(object)
#' print(x)
#'
#' ## cell_data_set ===
#' object <- cell_data_set
#' x <- cellCountsPerCluster(object)
#' print(x)
NULL



#' @rdname cellCountsPerCluster
#' @name cellCountsPerCluster
#' @importFrom bioverbs cellCountsPerCluster
#' @usage cellCountsPerCluster(object, ...)
#' @export
NULL



## nolint start
## dplyr approach:
## > x <- group_by(x, !!!syms(cols))
## > x <- summarize(x, n = n())
## > x <- ungroup(x)
## > x <- arrange(x, !!!syms(cols))
## > x <- group_by(x, !!sym("ident"))
## > x <- mutate(x, ratio = !!sym("n") / sum(!!sym("n")))
## nolint end



## Updated 2019-09-01.
`cellCountsPerCluster,SingleCellExperiment` <-  # nolint
    function(object, interestingGroups = NULL) {
        validObject(object)
        assert(.hasClusters(object))
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        x <- metrics(object, return = "DataFrame")
        ## Generate a cluster/sample contingency table from the metrics.
        tbl <- table(x[["ident"]], x[["sampleName"]])


        cols <- unique(c(
            "ident", "sampleName",
            interestingGroups, "interestingGroups"
        ))
        assert(isSubset(cols, colnames(x)))
        x <- x[, cols, drop = FALSE]
        f <- .group(x)
        split <- split(x = x, f = f)
        n <- vapply(
            X = split,
            FUN = nrow,
            FUN.VALUE = integer(1L),
            USE.NAMES = FALSE
        )
        ## Sample-level data.
        sl <- unique(x)
        rownames(sl) <- NULL
        sl <- sl
        ## Arrange by cluster then sample.
        sl <- sl[order(sl[["ident"]], sl[["sampleName"]]), , drop = FALSE]

        ## FIXME Rethink this approach...

        ## Generate the contingency table.
        tbl <- table(sl[["ident"]], sl[["sampleName"]])



        x <- arrange(x, !!!syms(cols))
        x <- group_by(x, !!sym("ident"))
        x <- mutate(x, ratio = !!sym("n") / sum(!!sym("n")))
        as(x, "DataFrame")
    }



#' @rdname cellCountsPerCluster
#' @export
setMethod(
    f = "cellCountsPerCluster",
    signature = signature("SingleCellExperiment"),
    definition = `cellCountsPerCluster,SingleCellExperiment`
)



## Updated 2019-07-31.
`cellCountsPerCluster,Seurat` <-  # nolint
    `cellCountsPerCluster,SingleCellExperiment`



#' @rdname cellCountsPerCluster
#' @export
setMethod(
    f = "cellCountsPerCluster",
    signature = signature("Seurat"),
    definition = `cellCountsPerCluster,Seurat`
)



## Updated 2019-08-02.
`cellCountsPerCluster,cell_data_set` <-  # nolint
    `cellCountsPerCluster,SingleCellExperiment`



#' @rdname cellCountsPerCluster
#' @export
setMethod(
    f = "cellCountsPerCluster",
    signature = signature("cell_data_set"),
    definition = `cellCountsPerCluster,cell_data_set`
)
