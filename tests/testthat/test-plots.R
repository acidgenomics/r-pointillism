context("Plot functions")

sce <- as(seurat, "SingleCellExperiment")
objects <- list(
    SingleCellExperiment = sce,
    seurat = seurat
)
genes <- head(rownames(seurat))



# plotCellTypesPerCluster ======================================================
with_parameters_test_that(
    "plotCellTypesPerCluster", {
        invisible(capture.output(
            list <- plotCellTypesPerCluster(
                object = object,
                markers = known_markers_small
            )
        ))
        expect_type(list, "list")
        expect_s3_class(list[[1L]][[1L]], "ggplot")
    },
    object = objects
)



# plotFeature ==================================================================
with_parameters_test_that(
    "plotFeature", {
        p <- plotFeature(object, features = c("PC1", "PC2"))
        expect_s3_class(p, "ggplot")
    },
    object = objects
)



# plotCounts ===================================================================
with_parameters_test_that(
    "plotCounts", {
        # Dot.
        p <- plotCounts(
            object = object,
            genes = genes,
            geom = "dot"
        )
        expect_s3_class(p, "ggplot")

        # Violin.
        p <- plotCounts(
            object = object,
            genes = genes,
            geom = "violin"
        )
        expect_s3_class(p, "ggplot")
    },
    object = objects
)

with_parameters_test_that(
    "plotDot", {
        p <- plotDot(object, genes = genes)
        expect_s3_class(p, "ggplot")
    },
    object = objects
)

with_parameters_test_that(
    "plotViolin", {
        p <- plotViolin(object, genes = genes)
        expect_s3_class(p, "ggplot")
    },
    object = objects
)




# plotMarker ===================================================================
with_parameters_test_that(
    "plotMarker", {
        expression <- methodFormals("plotMarker", "seurat") %>%
            .[["expression"]] %>%
            as.character() %>%
            .[-1L]
        invisible(lapply(expression, function(expression) {
            p <- plotMarker(
                object = object,
                genes = genes,
                expression = expression
            )
            expect_s3_class(p, "ggplot")
        }))
    },
    object = objects
)

with_parameters_test_that(
    "plotKnownMarkers", {
        invisible(capture.output(
            p <- plotKnownMarkers(
                object = object,
                markers = known_markers_small
            )
        ))
        expect_type(p, "list")
    },
    object = objects
)

with_parameters_test_that(
    "plotTopMarkers", {
        markers <- head(all_markers_small, n = 2L)
        invisible(capture.output(
            x <- plotTopMarkers(
                object = object,
                markers = markers
            )
        ))
        expect_type(x, "list")
        expect_s3_class(x[[1L]][[1L]], "ggplot")
    },
    object = objects
)



# plotPCElbow ==================================================================
# We're testing here to ensure that our seurat matches pbmc_small.
with_parameters_test_that(
    "plotPCElbow", {
        x <- plotPCElbow(object)
        expect_identical(x, seq)
    },
    object = list(
        seurat,
        Seurat::pbmc_small
    ),
    seq = list(
        seq_len(11L),
        seq_len(11L)
    )
)



# plotReducedDim ===============================================================
with_parameters_test_that(
    "plotReducedDim", {
        p <- plotReducedDim(
            object = object,
            reducedDim = "TSNE",
            pointsAsNumbers = TRUE,
            dark = TRUE,
            label = FALSE
        )
        expect_s3_class(p, "ggplot")
    },
    object = objects
)

# PCA
with_parameters_test_that(
    "plotPCA", {
        p <- plotPCA(object)
        expect_s3_class(p, "ggplot")
    },
    object = objects
)

# t-SNE
with_parameters_test_that(
    "plotTSNE", {
        p <- plotTSNE(object)
        expect_s3_class(p, "ggplot")
    },
    object = objects
)

# UMAP
with_parameters_test_that(
    "plotUMAP", {
        p <- plotTSNE(object)
        expect_s3_class(p, "ggplot")
    },
    object = objects
)
