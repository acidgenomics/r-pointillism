## Updated 2019-07-31.
.hasClusters <- function(object) {
    tryCatch(
        expr = is.factor(clusters(object)),
        error = function(e) FALSE
    )
}



## Updated 2019-07-31.
.hasDesignFormula <- function(object) {
    all(
        is(object, "SingleCellExperiment"),
        is.factor(object[["group"]]),
        is.matrix(metadata(object)[["design"]])
    )
}



## Updated 2019-07-31.
.hasMultipleSamples <- function(object) {
    length(sampleNames(object)) > 1L
}



## Consider moving this to goalie package.
## Updated 2019-07-31.
.isBPPARAM <- function(object) {
    all(
        identical(
            attributes(class(object))[["package"]],
            "BiocParallel"
        ),
        grepl("Param$", class(object))
    )
}
