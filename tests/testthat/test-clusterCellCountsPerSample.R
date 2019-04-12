context("clusterCellCountsPerSample")

with_parameters_test_that(
    "clusterCellCountsPerSample", {
        x <- clusterCellCountsPerSample(object)
        expect_is(x, "grouped_df")
        expect_identical(group_vars(x), "sampleName")
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
        SingleCellExperiment = sce,
        seurat = seurat_small
    )
)
