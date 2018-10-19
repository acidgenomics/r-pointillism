context("Differential Expression Functions")

data(seurat_small, envir = environment())

# Compare expression in cluster 3 relative to 2.
ident <- clusterID(seurat_small)
numerator <- names(ident)[ident == "3"]
denominator <- names(ident)[ident == "2"]

stopifnot(length(intersect(numerator, colnames(object))) > 0L)
stopifnot(length(intersect(denominator, colnames(object))) > 0L)

# Coerce back to seurat to match.
seurat_small <- as(sce_small, "seurat")



# diffExp ======================================================================
with_parameters_test_that(
    "diffExp", {
        # edgeR.
        x <- diffExp(
            object = object,
            numerator = numerator,
            denominator = denominator,
            caller = "edgeR"
        )
        expect_s4_class(x, "DGELRT")

        # DESeq2. Slow for large datasets.
        x <- diffExp(
            object = object,
            numerator = numerator,
            denominator = denominator,
            caller = "DESeq2"
        )
        expect_s4_class(x, "DESeqResults")
    },
    object = list(
        SingleCellExperiment = sce_small,
        seurat = seurat_small
    )
)



# findMarkers ==================================================================
with_parameters_test_that(
    "findMarkers", {
        # edgeR.
        x <- findMarkers(object, caller = "edgeR")
        expect_is(x, "list")
        invisible(lapply(
            X = x,
            FUN = function(x) {
                expect_is(x, "DGELRT")
            }
        ))

        # DESeq2.
        x <- findMarkers(object, caller = "DESeq2")
        expect_is(x, "list")
        invisible(lapply(
            X = x,
            FUN = function(x) {
                expect_is(x, "DESeqResults")
            }
        ))
    },
    object = list(
        SingleCellExperiment = sce_small
        # seurat = seurat_small
    )
)



# runZinbwave ==================================================================
test_that("runZinbwave", {
    # edgeR
    x <- runZinbwave(
        Y = sce_small,
        caller = "edgeR",
        recalculate = TRUE
    )

    # DESeq2
    x <- runZinbwave(
        Y = sce_small,
        caller = "DESeq2",
        recalculate = TRUE
    )
})
