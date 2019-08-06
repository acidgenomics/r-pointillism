#' @name metrics
#' @inherit basejump::metrics
#' @note Updated 2019-08-02.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(
#'     cell_data_set,
#'     package = "acidtest"
#' )
#'
#' ## cell_data_set ====
#' object <- cell_data_set
#' x <- metrics(object)
#' print(x)
#'
#' x <- metricsPerSample(object, fun = "mean")
#' print(x)
NULL



#' @rdname metrics
#' @name metrics
#' @importFrom bioverbs metrics
#' @usage metrics(object, ...)
#' @export
NULL



## This method will automatically add "ident" column and strip "cell" column.
## Updated 2019-08-06.
`metrics,cell_data_set` <-  # nolint
    function(object, return) {
        validObject(object)
        ## Strip invalid columns from column data.
        colData(object)[c("cell", "ident", "Size_Factor")] <- NULL
        data <- metrics(
            object = as(object, "SingleCellExperiment"),
            return = return
        )
        ident <- tryCatch(
            expr = clusters(object),
            error = function(e) NULL
        )
        if (is.factor(ident)) {
            data[["ident"]] <- ident
        }
        data
    }

f <- methodFormals(
    f = "metrics",
    signature = "SingleCellExperiment",
    package = "basejump"
)
formals(`metrics,cell_data_set`)[["return"]] <- f[["return"]]



#' @rdname metrics
#' @export
setMethod(
    f = "metrics",
    signature = signature("cell_data_set"),
    definition = `metrics,cell_data_set`
)
