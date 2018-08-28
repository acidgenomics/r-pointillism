.assertHasDesignFormula <- function(object) {
    stopifnot(is(object, "SingleCellExperiment"))
    assert_is_factor(object[["group"]])
    assert_is_matrix(metadata(object)[["design"]])
}



.assertHasIdent <- function(object) {
    assert_is_subset("ident", colnames(colData(object)))
}



.assertHasZinbwave <- function(object) {
    stopifnot(.hasZinbwave(object))
}



.assertIsKnownMarkers <- function(object) {
    # Require a tibble.
    assert_is_tbl_df(object)
    # Require grouping by `cellType` column.
    # Can use `attr()` instead of `group_vars()` here.
    assert_are_identical(group_vars(object), "cellType")
    # Require that there are genes.
    assert_has_rows(object)
}



.assertIsKnownMarkersDetected <- function(object) {
    .assertIsKnownMarkers(object)
    requiredCols <- c(
        "cellType",  # bcbio
        "cluster",   # Seurat
        "geneID",    # bcbio
        "geneName",  # bcbio
        "avgLogFC",  # Seurat v2.1
        "padj"       # Seurat v2.1
    )
    assert_is_subset(requiredCols, colnames(object))
}



.assertIsSanitizedMarkers <- function(object) {
    stopifnot(.isSanitizedMarkers(object))
}



.hasZinbwave <- function(object) {
    stopifnot(is(object, "SingleCellExperiment"))
    # Require `counts` to always be slotted.
    stopifnot("counts" %in% assayNames(object))
    all(c("normalizedValues", "weights") %in% assayNames(object))
}



.isSanitizedMarkers <- function(
    object,
    package = "Seurat"
) {
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
        # Check for `Seurat::FindAllMarkers()` return.
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
