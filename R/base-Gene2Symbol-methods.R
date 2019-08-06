#' @name Gene2Symbol
#' @inherit basejump::Gene2Symbol
#' @keywords internal
#' @note Updated 2019-08-02.
#'
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
#' x <- Gene2Symbol(object)
#' head(x)
#' table(x)
#'
#' ## cell_data_set ====
#' object <- cell_data_set
#' x <- Gene2Symbol(object)
#' head(x)
#' table(x)
NULL



#' @rdname Gene2Symbol
#' @name Gene2Symbol
#' @importFrom basejump Gene2Symbol
#' @usage Gene2Symbol(object, ...)
#' @export
NULL



## Updated 2019-08-02.
`Gene2Symbol,Seurat` <-  # nolint
    function(object, ...) {
        Gene2Symbol(
            object = as(object, "SummarizedExperiment"),
            ...
        )
    }



#' @rdname Gene2Symbol
#' @export
setMethod(
    f = "Gene2Symbol",
    signature = signature("Seurat"),
    definition = `Gene2Symbol,Seurat`
)



## Updated 2019-08-02.
`Gene2Symbol,cell_data_set` <-  # nolint
    function(object, ...) {
        assert(isSubset("gene_short_name", colnames(rowData(object))))
        df <- DataFrame(
            geneID = rownames(object),
            geneName = rowData(object)[["gene_short_name"]],
            row.names = rownames(object)
        )
        Gene2Symbol(object = df, ...)
    }



#' @rdname Gene2Symbol
#' @export
setMethod(
    f = "Gene2Symbol",
    signature = signature("cell_data_set"),
    definition = `Gene2Symbol,cell_data_set`
)
