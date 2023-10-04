## Updated 2019-08-23.
.getSeuratStash <- function(object, name) {
    assert(
        is(object, "Seurat"),
        isString(name)
    )
    misc <- slot(object, name = "misc")
    ## Early return if the `misc` slot is `NULL`.
    if (is.null(misc)) {
        return(NULL)
    }
    ## Look first directly in `object@misc` slot.
    x <- misc[[name]]
    if (!is.null(x)) {
        return(x)
    }
    ## Next, handle legacy `bcbio` stash list inside `object@misc`.
    ## As of v0.1.3, stashing directly into `object@misc`.
    if ("bcbio" %in% names(misc)) {
        x <- misc[["bcbio"]][[name]]
        if (!is.null(x)) {
            return(x)
        }
    }
    NULL
}



## Updated 2020-02-21.
.seuratCommand <-
    function(object, name, assay = NULL) {
        assert(
            is(object, "Seurat"),
            isString(name),
            isString(assay, nullOk = TRUE)
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
#' @note Updated 2021-09-03.
#' @noRd
#'
#' @seealso
#' - `Seurat::Idents()`
#' - `Seurat:::Idents.Seurat`
.seuratWhichIdents <- function(object) {
    data <- slot(object, name = "meta.data")
    idents <- Idents(object)
    assert(
        is.data.frame(data),
        is.factor(idents),
        identical(rownames(data), names(idents))
    )
    keep <- grepl("res\\.[.0-9]+$", colnames(data))
    assert(
        any(keep),
        msg = sprintf(
            "Failed to detect any resolutions in '%s'.",
            "object@meta.data"
        )
    )
    data <- data[, keep, drop = FALSE]
    ## Now check against cluster idents currently returned by `Idents()`, which
    ## are internally stashed in `object@active.ident`.
    keep <- bapply(
        X = data,
        y = idents,
        FUN = function(x, y) {
            assert(
                is.factor(x),
                is.factor(y)
            )
            x <- unname(x)
            x <- droplevels(x)
            y <- unname(y)
            y <- droplevels(y)
            identical(x, y)
        }
    )
    assert(
        any(keep),
        msg = sprintf(
            "Failed to match '%s' in '%s'.",
            "Idents()", "object@meta.data"
        )
    )
    col <- names(keep)[keep]
    ## Inform the user about multiple identical resolutions, which can happen
    ## with low complexity samples.
    if (!isString(col)) {
        alertWarning(sprintf(
            "Multiple resolutions matched: %s",
            toInlineString(col, n = 5L)
        ))
        col <- col[[1L]]
    }
    col
}
