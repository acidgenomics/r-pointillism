context("Cluster Statistics Functions")

# basejump SCE example doesn't contain clustering data.
# data(sce_small, package = "basejump", envir = environment())

data(seurat_small, known_markers_small, envir = environment())
sce_small <- as(seurat_small, "SingleCellExperiment")



# cellCountsPerCluster ========================================================
with_parameters_test_that(
    "cellCountsPerCluster", {
        x <- cellCountsPerCluster(object)
        expect_is(x, "grouped_df")
        expect_identical(dplyr::group_vars(x), "ident")
    },
    object = list(
        SingleCellExperiment = sce_small,
        seurat = seurat_small
    )
)



# cellTypesPerCluster ==========================================================
with_parameters_test_that(
    "cellTypesPerCluster", {
        x <- cellTypesPerCluster(object)
        expect_is(x, "grouped_df")
        expect_identical(dplyr::group_vars(x), "cluster")
        expect_identical(
            lapply(x, class),
            list(
                cluster = "factor",
                cellType = "factor",
                n = "integer",
                name = "character",
                geneID = "character",
                geneName = "character"
            )
        )
    },
    object = list(
        KnownSeuratMarkers = known_markers_small
    )
)



# clusterCellCountsPerSample ==================================================
with_parameters_test_that(
    "clusterCellCountsPerSample", {
        x <- clusterCellCountsPerSample(object)
        expect_is(x, "grouped_df")
        expect_identical(dplyr::group_vars(x), "sampleName")
        expect_identical(
            lapply(x, class),
            list(
                sampleName = "factor",
                ident = "factor",
                n = "integer",
                ratio = "numeric"
            )
        )
    },
    object = list(
        SingleCellExperiment = sce_small,
        seurat = seurat_small
    )
)
