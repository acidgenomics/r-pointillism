#' Cell-Cycle Markers
#'
#' Currently designed for internal use by the pointillism package.
#'
#' @family Marker Generators
#' @include markers-internal.R
#' @export
#'
#' @inheritParams CellTypeMarkers
#'
#' @return `CellCycleMarkers`.
#' @examples
#' \dontrun{
#' # Google Sheets method (requires OAuth).
#' x <- CellCycleMarkers(
#'     gs = "1qA5ktYeimNGpZF1UPSQZATbpzEqgyxN6daoMOjv6YYw",
#'     ws = "Homo_sapiens",
#'     organism = "Homo sapiens",
#'     ensemblRelease = 92L
#' )
#' }
CellCycleMarkers <-  # nolint
    function() {
        do.call(
            what = .getMarkers,
            args = matchArgsToDoCall(args = list(class = "CellCycleMarkers"))
        )
    }
f <- formals(.getMarkers)
f <- f[setdiff(names(f), "class")]
formals(CellCycleMarkers) <- f
