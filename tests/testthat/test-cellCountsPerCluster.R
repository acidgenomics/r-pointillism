context("cellCountsPerCluster")

with_parameters_test_that(
    "cellCountsPerCluster", {
        x <- cellCountsPerCluster(object)
        expect_is(x, "grouped_df")
        expect_identical(group_vars(x), "ident")
    },
    object = list(
        SingleCellExperiment = sce,
        seurat = seurat
    )
)
