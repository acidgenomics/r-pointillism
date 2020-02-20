#' Run Seurat
#'
#' @export
#' @note Updated 2020-02-20.
#'
#' @section reticulate:
#'
#' Check `Sys.getenv("WORKON_HOME")` path, which is the current approach used
#' in `Rprofile.site`.
#'
#' @seealso
#' - https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
#' - https://github.com/satijalab/seurat/wiki
runSeurat <- function(
    object,
    virtualenv = "r-reticulate",
    regressCellCycle = TRUE
) {
    assert(
        requireNamespace("reticulate", quietly = TRUE),
        requireNamespace("Seurat", quietly = TRUE),
        is(object, "Seurat"),
        isString(virtualenv),
        isFlag(regressCellCycle)
    )
    reticulate::use_virtualenv(
        virtualenv = virtualenv,
        required = TRUE
    )
    assert(py_module_available(module = "umap"))

    if (isTRUE(regressCellCycle)) {
        ## Get CellCycleMarkers object required for cell-cycle regerssion.
        organism <- organism(object)
        data(
            cell_cycle_markers_list,
            package = "pointillism",
            envir = environment()
        )
        ccm <- cell_cycle_markers_list[[camelCase(organism)]]
        assert(is(ccm, "CellCycleMarkers"))
        g2mGenes <- as.character(ccm[["g2m"]][["geneName"]])
        sGenes <- as.character(ccm[["s"]][["geneName"]])
    }

    object <- Seurat::NormalizeData(object)
    object <- Seurat::FindVariableFeatures(object)
    object <- Seurat::ScaleData(object)

    object
}
