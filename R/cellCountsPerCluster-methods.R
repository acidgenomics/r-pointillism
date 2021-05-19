#' @name cellCountsPerCluster
#' @inherit AcidGenerics::cellCountsPerCluster
#' @note Updated 2021-03-03.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @return `DataFrame`.
#'
#' @examples
#' data(Seurat, package = "AcidTest")
#'
#' ## Seurat ====
#' object <- Seurat
#' x <- cellCountsPerCluster(object)
#' print(x)
NULL



## Updated 2019-10-30.
`cellCountsPerCluster,SingleCellExperiment` <-  # nolint
    function(object) {
        validObject(object)
        assert(.hasClusters(object))
        interestingGroups <- interestingGroups(object)
        x <- metrics(object, return = "DataFrame")
        ## Contingency table.
        tbl <- table(x[["ident"]], x[["sampleId"]])
        ## Get the number of cells per ident.
        nPerIdent <- rowSums(tbl)
        nPerIdent <- DataFrame(
            ident = names(nPerIdent),
            nPerIdent = as.integer(nPerIdent)
        )
        ## Summarize to cluster/sample-level metadata.
        cols <- unique(c(
            "ident", "sampleId", "sampleName",
            interestingGroups, "interestingGroups"
        ))
        x <- unique(x[, cols, drop = FALSE])
        rownames(x) <- NULL
        x <- x[order(x[["ident"]], x[["sampleId"]]), , drop = FALSE]
        ## Melt the contingency table into long format.
        melt <- melt(tbl, colnames = c("ident", "sampleId", "n"))
        ## Join and calculate ratio.
        x <- leftJoin(x, melt, by = c("ident", "sampleId"))
        x <- leftJoin(x, nPerIdent, by = "ident")
        x[["ratio"]] <- x[["n"]] / x[["nPerIdent"]]
        x
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
