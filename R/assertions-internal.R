.hasDesignFormula <- function(object) {
    all(
        is(object, "SingleCellExperiment"),
        is.factor(object[["group"]]),
        is.matrix(metadata(object)[["design"]])
    )
}



.hasIdent <- function(object) {
    "ident" %in% colnames(colData(object))
}



.hasMultipleSamples <- function(object) {
    length(sampleNames(object)) > 1L
}



.isBPPARAM <- function(object) {
    all(
        identical(
            attributes(class(object))[["package"]],
            "BiocParallel"
        ),
        grepl("Param$", class(object))
    )
}



.isKnownMarkers <- function(object) {
    all(
        # Require a tibble.
        is(object, "tbl_df"),
        # Require grouping by `cellType` column.
        # Can use `attr` instead of `group_vars` here.
        identical(group_vars(object), "cellType"),
        # Require that there are genes.
        hasRows(object)
    )
}



.isKnownMarkersDetected <- function(object) {
    .isKnownMarkers(object)
    requiredCols <- c(
        "cellType",  # bcbio
        "cluster",   # Seurat
        "geneID",    # bcbio
        "geneName",  # bcbio
        "avgLogFC",  # Seurat v2.1
        "padj"       # Seurat v2.1
    )
    isSubset(requiredCols, colnames(object))
}



.isSanitizedMarkers <- function(object) {
    assert(.isSanitizedMarkers(object))
}



.isSanitizedMarkers <- function(object, package = "Seurat") {
    package <- match.arg(package)

    # General checks -----------------------------------------------------------
    if (!is(object, "grouped_df")) {
        return(FALSE)
    } else if (
        is.null(attr(object, "vars")) ||
        attr(object, "vars") != "cluster"
    ) {
        return(FALSE)
    } else if (!"geneID" %in% colnames(object)) {
        return(FALSE)
    }

    # Package-specific checks --------------------------------------------------
    if (package == "Seurat") {
        # Check for `Seurat::FindAllMarkers` return.
        # These columns are output in an inconsistent format, so we'll sanitize
        # into lowerCamelCase.
        seuratBlacklist <- c(
            "avg_diff",   # Legacy, now "avg_logFC"
            "avg_logFC",  # Renamed in v2.1
            "gene",
            "p_val",      # We'll rename to pvalue, matching DESeq2
            "p_val_adj",  # New in v2.1, we'll rename to padj, matching DESeq2
            "pct.1",
            "pct.2"
        )
        if (any(seuratBlacklist %in% colnames(object))) {
            return(FALSE)
        } else {
            return(TRUE)
        }
    }
}
