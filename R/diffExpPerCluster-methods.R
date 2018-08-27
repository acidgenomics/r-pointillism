#' Differential Expression per Cluster
#'
#' @note Cluster identity (`ident`) must be defined in `colData()` for this
#'   function to work.
#'
#' @name diffExpPerCluster
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param ... Passthrough arguments to [diffExp()].
#'
#' @return `list` containing:
#' - `caller = "edgeR"`: `DGELRT`.
#' - `caller = "DESeq2"`: `DESeqResults`.
#'
#' @examples
#' x <- suppressMessages(diffExpPerCluster(
#'     object = sce_small,
#'     group = "group",
#'     numerator = "group2",
#'     denominator = "group1",
#'     caller = "edgeR"
#' ))
#' class(x)
#' lapply(x, class)
NULL



#' @rdname diffExpPerCluster
#' @export
setMethod(
    "diffExpPerCluster",
    signature("SingleCellExperiment"),
    function(
        object,
        group,
        numerator,
        denominator,
        ...
    ) {
        # Object must contain pre-calculate ZINB weights.
        .assertHasZinbwave(object)
        # group
        assert_is_a_string(group)
        assert_is_subset(group, colnames(colData(object)))
        groupings <- colData(object)[[group]]
        assert_is_factor(groupings)
        # numerator
        assert_is_a_string(numerator)
        assert_is_subset(numerator, levels(groupings))
        # denominator
        assert_is_a_string(denominator)
        assert_is_subset(denominator, levels(groupings))
        assert_are_disjoint_sets(numerator, denominator)
        # Get the cluster identities.
        ident <- clusterID(object)
        assert_is_factor(ident)
        clusters <- levels(ident)
        stopifnot(length(clusters) >= 2L)
        message(paste(length(clusters), "clusters detected"))
        # Loop across each cluster and perform pairwise DE based on the single
        # group of interest defined.
        # Consider adding a skip step here for clusters with very few cells.
        list <- lapply(
            X = clusters,
            FUN = function(cluster) {
                message(paste("Cluster", cluster, "===="))
                # Subset the cells by cluster.
                cells <- colnames(object)[which(ident == cluster)]
                assert_is_non_empty(cells)
                subset <- object[, cells]
                # Ensure that both the numerator and denominator are defined.
                groupings <- droplevels(colData(subset)[[group]])
                numerator <- colnames(subset)[which(groupings == numerator)]
                denominator <- colnames(subset)[which(groupings == denominator)]
                if (
                    length(numerator) == 0L ||
                    length(denominator) == 0L
                ) {
                    .warnBadContrast()
                    return(NULL)
                }
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
)
