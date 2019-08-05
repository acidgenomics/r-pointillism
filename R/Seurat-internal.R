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



.seuratCommandParam <-
    function(object, name, param, assay = NULL) {
        cmd <- .seuratCommand(object = object, name = name, assay = assay)
        assert(.hasSlot(cmd, "params"))
        cmd@params[[param]]
    }



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
