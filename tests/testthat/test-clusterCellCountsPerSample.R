## FIXME Redundant with `cellCountsPerCluster()`?

context("clusterCellCountsPerSample")

with_parameters_test_that(
    "clusterCellCountsPerSample", {
        colData(object)[["sampleID"]] <- as.factor(paste0("sample", seq(2L)))
        colData(object)[["sampleName"]] <- colData(object)[["sampleID"]]
        object <- clusterCellCountsPerSample(object)
        expect_s4_class(object, "DataFrame")
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
