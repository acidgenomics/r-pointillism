context("Differential Expression")

data(seurat_small, envir = environment())
sce <- as(seurat_small, "SingleCellExperiment")

# Compare expression in cluster 3 relative to 2.
object <- seurat_small
ident <- clusterID(object)
numerator <- names(ident)[ident == "3"]
denominator <- names(ident)[ident == "2"]
expect_true(length(intersect(numerator, colnames(object))) > 0L)
expect_true(length(intersect(denominator, colnames(object))) > 0L)
rm(object)



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
        # Expecting warning about degenerate design matrix.
        x <- suppressWarnings(
            diffExp(
                object = object,
                numerator = numerator,
                denominator = denominator,
                caller = "DESeq2"
            )
        )
        expect_s4_class(x, "DESeqResults")
    },
    object = list(
        SingleCellExperiment = sce,
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

        # DESeq2. Slow for large datasets.
        # Expecting warning about degenerate design matrix.
        x <- suppressWarnings(
            findMarkers(object, caller = "DESeq2")
        )
        expect_is(x, "list")
        invisible(lapply(
            X = x,
            FUN = function(x) {
                expect_is(x, "DESeqResults")
            }
        ))
    },
    object = list(
        SingleCellExperiment = sce,
        seurat = seurat_small
    )
)
