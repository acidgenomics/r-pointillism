#' Show an object
#' @name show
#' @author Michael Steinbuagh
#' @inherit methods::show title description details params
#' @note Updated 2020-02-20.
#'
#' @return No return.
#'
#' @examples
#' data(cell_cycle_markers_list, cell_type_markers_list)
#'
#' ## CellCycleMarkers ====
#' object <- cell_cycle_markers_list[[1L]]
#' show(object)
#'
#' ## CellTypeMarkers ====
#' object <- cell_type_markers_list[[1L]]
#' show(object)
NULL



## Updated 2020-02-20.
`show,CellCycleMarkers` <-  # nolint
    function(object) {
        validObject(object)
        ## Include the organism information.
        organism <- metadata(object)[["organism"]]
        release <- metadata(object)[["ensemblRelease"]]
        lengths <- nrow(object)
        genes <- vapply(
            X = object,
            FUN = function(x) {
                x <- as.character(x[["geneName"]])
                x <- sort(x)
                x <- c(head(x, n = 2L), "...", tail(x, n = 2L))
                paste(x, collapse = " ")
            },
            FUN.VALUE = character(1L),
            USE.NAMES = TRUE
        )
        return <- c(
            class(object),
            paste0(organism, " (Ensembl ", release, ")"),
            paste0(names(genes), "(", lengths, "): ", genes)
        )
        cat(return, sep = "\n")
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("CellCycleMarkers"),
    definition = `show,CellCycleMarkers`
)



## Updated 2020-02-20.
`show,CellTypeMarkers` <-  # nolint
    function(object) {
        validObject(object)

        ## Include the organism information.
        organism <- metadata(object)[["organism"]]
        release <- metadata(object)[["ensemblRelease"]]

        ## Include the gene lengths per phase.
        lengths <- nrow(object)
        genes <- vapply(
            X = object,
            FUN = function(x) {
                x <- as.character(x[["geneName"]])
                x <- sort(x)
                x <- c(head(x, n = 2L), "...", tail(x, n = 2L))
                paste(x, collapse = " ")
            },
            FUN.VALUE = character(1L),
            USE.NAMES = TRUE
        )

        return <- c(
            class(object),
            paste0(organism, " (Ensembl ", release, ")"),
            paste0(names(genes), "(", lengths, "): ", genes)
        )

        cat(return, sep = "\n")
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("CellTypeMarkers"),
    definition = `show,CellTypeMarkers`
)
