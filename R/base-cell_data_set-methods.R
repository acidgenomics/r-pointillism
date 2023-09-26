## > #' Extend base S4 methods for `cell_data_set` class
## > #'
## > #' @name base-cell_data_set
## > #' @keywords internal
## > #' @note Updated 2021-10-13.
## > #'
## > #' @inheritParams AcidRoxygen::params
## > #'
## > #' @return Varies, depending on the generic.
## > NULL



## > ## Updated 2019-08-02.
## > `GeneToSymbol,cell_data_set` <-  # nolint
## >     function(object, ...) {
## >         assert(isSubset("gene_short_name", colnames(rowData(object))))
## >         df <- DataFrame(
## >             "geneId" = rownames(object),
## >             "geneName" = rowData(object)[["gene_short_name"]],
## >             row.names = rownames(object)
## >         )
## >         GeneToSymbol(object = df, ...)
## >     }
## >
## > #' @rdname base-cell_data_set
## > #' @export
## > setMethod(
## >     f = "GeneToSymbol",
## >     signature = signature(object = "cell_data_set"),
## >     definition = `GeneToSymbol,cell_data_set`
## > )



## > ## Note that monocle3 currently defines the generic using "x" instead of
## > ## "object", and requires "reduction_method" in definition.
## > ##
## > ## For `reduction`, note that positional scalar works.
## > ##
## > ## Updated 2019-08-02.
## > `clusters,cell_data_set` <-  # nolint
## >     function(object, reduction) {
## >         validObject(object)
## >         assert(isScalar(reduction))
## >         monocle3::clusters(
## >             x = object,
## >             reduction_method = reduction
## >         )
## >     }
## >
## > f <- methodFormals(
## >     f = "clusters",
## >     signature = signature(x = "cell_data_set"),
## >     package = "monocle3"
## > )
## > formals(`clusters,cell_data_set`)[["reduction"]] <-
## >     f[["reduction_method"]]
## >
## > #' @rdname clusters
## > #' @export
## > setMethod(
## >     f = "clusters",
## >     signature = signature(object = "cell_data_set"),
## >     definition = `clusters,cell_data_set`
## > )



## > ## Updated 2020-01-30.
## > `logcounts,cell_data_set` <-  # nolint
## >     function(object) {
## >         monocle3::normalized_counts(
## >             cds = object,
## >             norm_method = "log",
## >             pseudocount = 1L
## >         )
## >     }
## >
## > #' @rdname base-cell_data_set
## > #' @export
## > setMethod(
## >     f = "logcounts",
## >     signature = signature(object = "cell_data_set"),
## >     definition = `logcounts,cell_data_set`
## > )



## > ## This method will automatically add "ident" and strip "cell" column.
## > ## Updated 2019-08-06.
## > `metrics,cell_data_set` <-  # nolint
## >     function(object, return) {
## >         validObject(object)
## >         ## Strip invalid columns from column data.
## >         colData(object)[c("cell", "ident", "Size_Factor")] <- NULL
## >         data <- metrics(
## >             object = as(object, "SingleCellExperiment"),
## >             return = return
## >         )
## >         ident <- tryCatch(
## >             expr = clusters(object),
## >             error = function(e) NULL
## >         )
## >         if (is.factor(ident)) {
## >             data[["ident"]] <- ident
## >         }
## >         data
## >     }
## >
## > f <- methodFormals(
## >     f = "metrics",
## >     signature = "SingleCellExperiment",
## >     package = "basejump"
## > )
## > formals(`metrics,cell_data_set`)[["return"]] <- f[["return"]]
## >
## > #' @rdname base-cell_data_set
## > #' @export
## > setMethod(
## >     f = "metrics",
## >     signature = signature(object = "cell_data_set"),
## >     definition = `metrics,cell_data_set`
## > )



## > ## Updated 2020-01-30.
## > `normalize,cell_data_set` <-  # nolint
## >     function(object) {
## >         alert(
## >             "Normalizing with {.pkg monocle3}::{.fun preprocess_cds}."
## >         )
## >         monocle3::preprocess_cds(
## >             cds = object,
## >             method = "PCA",
## >             norm_method = "log",
## >             pseudo_count = 1L,
## >             scaling = TRUE,
## >             verbose = TRUE
## >         )
## >     }
## >
## > #' @rdname normalize
## > #' @export
## > setMethod(
## >     f = "normalize",
## >     signature = signature(object = "cell_data_set"),
## >     definition = `normalize,cell_data_set`
## > )



## > ## Updated 2020-01-30.
## > `normcounts,cell_data_set` <-  # nolint
## >     function(object) {
## >         monocle3::normalized_counts(
## >             cds = object,
## >             norm_method = "size_only",
## >             pseudocount = NULL
## >         )
## >     }
## >
## > #' @rdname base-cell_data_set
## > #' @export
## > setMethod(
## >     f = "normcounts",
## >     signature = signature(object = "cell_data_set"),
## >     definition = `normcounts,cell_data_set`
## > )



## > ## Updated 2019-08-06.
## > `sizeFactors,cell_data_set` <-  # nolint
## >     function(object) {
## >         colData(object)[["Size_Factor"]]
## >     }
## >
## > #' @rdname base-cell_data_set
## > #' @export
## > setMethod(
## >     f = "sizeFactors",
## >     signature = signature(object = "cell_data_set"),
## >     definition = `sizeFactors,cell_data_set`
## > )



## > ## Updated 2019-08-06.
## > `sizeFactors<-,cell_data_set,ANY` <-  # nolint
## >     function(object, value) {
## >         if (!is.null(value)) {
## >             assert(
## >                 all(!is.na(value)),
## >                 all(is.finite(value)),
## >                 all(value > 0L)
## >             )
## >             value <- unname(value)
## >         }
## >         colData(object)[["Size_Factor"]] <- value
## >         object
## >     }
## >
## > #' @rdname base-cell_data_set
## > #' @export
## > setReplaceMethod(
## >     f = "sizeFactors",
## >     signature = signature(
## >         object = "cell_data_set",
## >         value = "ANY"
## >     ),
## >     definition = `sizeFactors<-,cell_data_set,ANY`
## > )
