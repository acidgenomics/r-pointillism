## Consider adding pseudobulk support here in a future update.



#' @name diffExpPerCluster
#' @inherit AcidGenerics::diffExpPerCluster
#' @note Updated 2020-01-30.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams diffExp
#' @param group `character(1)`.
#'   Group of interest for differential expression per cluster. Must be a
#'   `factor` column in [`colData()`][SummarizedExperiment::colData].
#' @param ... Passthrough arguments to [diffExp()].
#'
#' @note Cluster identity (`ident`) must be defined in
#'   [`colData()`][SummarizedExperiment::colData] for this function to work.
#'
#' @return `list` containing:
#' - `caller = "edgeR"`: `DGELRT`.
#' - `caller = "DESeq2"`: `DESeqResults`.
#'
#' @examples
#' data(Seurat, package = "AcidTest")
#'
#' ## Seurat ====
#' object <- Seurat
#' group <- factor(c("group1", "group2"))
#' colData(object)$group <- group
#' suppressMessages({
#'     x <- diffExpPerCluster(
#'         object = object,
#'         group = "group",
#'         numerator = "group2",
#'         denominator = "group1",
#'         caller = "edgeR"
#'     )
#' })
#' class(x)
#' lapply(x, class)
NULL



## Updated 2020-01-30.
`diffExpPerCluster,SingleCellExperiment` <-  # nolint
    function(
        object,
        group,
        numerator,
        denominator,
        ...
    ) {
        h1("{.fun diffExpPerCluster}")
        object <- as(object, "SingleCellExperiment")
        ## group
        assert(
            isString(group),
            isSubset(group, colnames(colData(object)))
        )
        groupings <- colData(object)[[group]]
        assert(is.factor(groupings))
        ## numerator
        assert(
            isString(numerator),
            isSubset(numerator, levels(groupings))
        )
        ## denominator
        assert(
            isString(denominator),
            isSubset(denominator, levels(groupings)),
            areDisjointSets(numerator, denominator)
        )
        ## Get the cluster identities.
        ident <- clusters(object)
        assert(is.factor(ident))
        clusters <- levels(ident)
        assert(length(clusters) >= 2L)
        alertInfo(sprintf("%d clusters detected.", length(clusters)))
        ## Loop across each cluster and perform pairwise DE based on the single
        ## group of interest defined.
        ## Consider adding a skip step here for clusters with very few cells.
        list <- lapply(
            X = clusters,
            FUN = function(cluster) {
                h2(sprintf("Cluster %s", as.character(cluster)))
                ## Subset the cells by cluster.
                cells <- colnames(object)[which(ident == cluster)]
                assert(hasLength(cells))
                subset <- object[, cells]
                ## Ensure that both the numerator and denominator are defined.
                groupings <- droplevels(colData(subset)[[group]])
                numerator <- colnames(subset)[which(groupings == numerator)]
                denominator <- colnames(subset)[which(groupings == denominator)]
                diffExp(
                    object = object,
                    numerator = numerator,
                    denominator = denominator,
                    ...
                )
            }
        )
        names(list) <- clusters
        list
    }



#' @rdname diffExpPerCluster
#' @export
setMethod(
    f = "diffExpPerCluster",
    signature = signature("SingleCellExperiment"),
    definition = `diffExpPerCluster,SingleCellExperiment`
)



## Updated 2019-07-31.
`diffExpPerCluster,Seurat` <-  # nolint
    `diffExpPerCluster,SingleCellExperiment`



#' @rdname diffExpPerCluster
#' @export
setMethod(
    f = "diffExpPerCluster",
    signature = signature("Seurat"),
    definition = `diffExpPerCluster,Seurat`
)
