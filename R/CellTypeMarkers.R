#' Cell-Type Markers
#'
#' Local spreadsheets (CSV, Excel) and Google Sheets are currently supported.
#'
#' Must contain `phase` and `geneID` columns.
#'
#' @family Marker Generators
#' @include markers-internal.R
#' @export
#'
#' @inheritParams general
#' @inheritParams googlesheets::gs_read
#' @param file `string` or `missing`. Gene markers file path (CSV or Excel).
#' @param gs `string` or `missing`. Google Sheets `sheet_key` identifier, which
#'   can be located with the [googlesheets::gs_ls()] function.
#' @param gene2symbol `gene2symbol` or `NULL`. Use the stashed gene-to-symbol
#'   mappings from a `SingleCellExperiment` or `seurat` object.
#'
#' @seealso
#' - [googlesheets::gs_key()].
#' - [googlesheets::gs_read()].
#'
#' @return `CellTypeMarkers`.
#' @examples
#' # Google Sheets method (requires OAuth).
#' \dontrun{
#' x <- CellTypeMarkers(
#'     gs = "1vGNU2CCxpaoTCLvzOxK1hf5gjULrf2-CpgCp9bOfGJ0",
#'     ws = "Homo_sapiens",
#'     organism = "Homo sapiens",
#'     ensemblRelease = 92L
#' )
#' }
CellTypeMarkers <-  # nolint
    function() {
        do.call(
            what = .getMarkers,
            args = matchArgsToDoCall(args = list(class = "CellTypeMarkers"))
        )
    }
f <- formals(.getMarkers)
f <- f[setdiff(names(f), "class")]
formals(CellTypeMarkers) <- f
