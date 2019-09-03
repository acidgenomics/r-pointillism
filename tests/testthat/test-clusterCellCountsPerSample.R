context("clusterCellCountsPerSample")

with_parameters_test_that(
    "clusterCellCountsPerSample", {
        x <- clusterCellCountsPerSample(object)
        expect_s4_class(x, "DataFrame")
        expect_identical(
            lapply(x, class),
            list(
                ident = "factor",
                sampleName = "factor",
                n = "integer",
                nPerIdent = "integer",
                ratio = "numeric"
            )
        )
    },
    object = list(
        SingleCellExperiment = sce,
        Seurat = seurat
    )
)
