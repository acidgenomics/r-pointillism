#' Import markers from Google Sheets
#'
#' Currently designed for internal use by the pointillism package.
#'
#' @name importMarkers
#'
#' @inheritParams basejump::params
#' @inheritParams basejump::makeGRanges
#' @inheritParams googlesheets::gs_read
#'
#' @param gs `character(1)`.
#'   Google Sheets `sheet_key` identifier, which can be located with the
#'   [googlesheets::gs_ls()] function.
#' @param gene2symbol `Gene2Symbol`.
#'   Gene-to-symbol mappings.
#'
#' @return `CellCycleMarkers` or `CellTypeMarkers`.
#'
#' @note Google Sheets requires OAuth authentication.
#'
#' @seealso
#' - `googlesheets::gs_key()`.
#' - `googlesheets::gs_read()`.
#'
#' @examples
#' ## Must be interactive.
#' if (isTRUE(interactive())) {
#'     googlesheets::gs_ls()
#'
#'     ## Gene-to-symbol mappings.
#'     gene2symbol <- makeGene2SymbolFromEnsembl("Homo sapiens")
#'
#'     ## Cell-cycle markers.
#'     x <- importCellCycleMarkersFromGoogle(
#'         gs = "1qA5ktYeimNGpZF1UPSQZATbpzEqgyxN6daoMOjv6YYw",
#'         ws = "Homo_sapiens",
#'         gene2symbol = gene2symbol
#'     )
#'
#'     ## Cell-type markers.
#'     x <- importCellTypeMarkersFromGoogle(
#'         gs = "1vGNU2CCxpaoTCLvzOxK1hf5gjULrf2-CpgCp9bOfGJ0",
#'         ws = "Homo_sapiens",
#'         gene2symbol = gene2symbol
#'     )
#' }
NULL



# Import data from Google Sheets.
.importFromGoogle <- function(gs, ws) {
    assert(
        isString(gs),
        isScalar(ws)
    )
    ss <- gs_key(gs)
    data <- gs_read(ss = ss, ws = ws)
    as(data, "DataFrame")
}



#' @describeIn importMarkers Import Cell-Cycle Markers from Google
#' @export
importCellCycleMarkersFromGoogle <-
    function(gs, ws, gene2symbol) {
        data <- .importFromGoogle(gs = gs, ws = ws)
        CellCycleMarkers(object = data, gene2symbol = gene2symbol)
    }



#' @describeIn importMarkers Cell-type markers.
#' @export
importCellTypeMarkersFromGoogle <-
    function(gs, ws, gene2symbol) {
        data <- .importFromGoogle(gs = gs, ws = ws)
        CellTypeMarkers(object = data, gene2symbol = gene2symbol)
    }
