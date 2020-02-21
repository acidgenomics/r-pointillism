## Updated 2020-02-21.
.seuratCommand <-
    function(object, name, assay = NULL) {
        assert(
            is(object, "Seurat"),
            isString(name),
            isString(assay, nullOK = TRUE)
        )
        if (is.null(assay)) {
            assay <- DefaultAssay(object)
            assert(isString(assay))
        }
        commands <- slot(object, "commands")
        assert(is.list(commands))
        cmd <- commands[[paste0(name, ".", assay)]]
        assert(is(cmd, "SeuratCommand"))
        cmd
    }



## Updated 2020-02-21.
.seuratCommandParam <-
    function(object, name, param, assay = NULL) {
        cmd <- .seuratCommand(object = object, name = name, assay = assay)
        assert(.hasSlot(cmd, "params"))
        cmd@params[[param]]
    }



## Updated 2020-02-21.
.seuratNormalizationMethod <- function(object, assay = NULL) {
    x <- .seuratCommandParam(
        object = object,
        name = "NormalizeData",
        param = "normalization.method",
        assay = assay
    )
    assert(isString(x))
    x
}



## Updated 2020-02-21.
.seuratScaleFactor <- function(object, assay = NULL) {
    x <- .seuratCommandParam(
        object = object,
        name = "NormalizeData",
        param = "scale.factor",
        assay = assay
    )
    assert(isNumber(x))
    x
}



#' Determine which cluster resolution maps to active cell idents
#'
#' Internally this checks `Idents()` return against internal metadata columns
#' slotted in `object@meta.data`.
#'
#' @note Updated 2020-02-21.
#' @noRd
#'
#' @seealso
#' - `Seurat::Idents()`
#' - `Seurat:::Idents.Seurat`
.seuratWhichIdents <- function(object) {
    data <- object@meta.data
    keep <- grepl("res\\.[.0-9]+$", colnames(data))
    if (!any(keep)) {
        stop("Failed to detect any resolutions in `object@meta.data`.")
    }
    data <- data[, keep, drop = FALSE]
    ## Now check against cluster idents currently returned by `Idents()`, which
    ## are internally stashed in `object@active.ident`.
    keep <- bapply(
        X = data,
        y = Idents(object),
        FUN = function(x, y) {
            identical(unname(x), unname(y))
        }
    )
    if (!any(keep)) {
        stop("Failed to match `Idents()` in `object@meta.data`.")
    }
    col <- names(keep)[keep]
    assert(isString(col))
    col
}
