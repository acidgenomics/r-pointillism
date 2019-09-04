## FIXME How is this different from `clusterCellCountsPerSample()`?

context("cellCountsPerCluster")

with_parameters_test_that(
    "cellCountsPerCluster", {
        ## Simulate multiple samples.
        samples <- factor(paste0("sample", seq_len(2L)))
        colData(object)[["sampleID"]] <- samples
        colData(object)[["sampleName"]] <- samples
        object <- cellCountsPerCluster(object)
        expect_s4_class(object, "DataFrame")
        expected <- DataFrame(
            ident = as.factor(rep(c("0", "1", "2"), each = 2L)),
            sampleID = samples,
            sampleName = samples,
            interestingGroups = samples,
            n = c(19L, 17L, 11L, 14L, 10L, 9L),
            nPerIdent = rep(c(36L, 25L, 19L), each = 2L)
        )
        expected[["ratio"]] <- expected[["n"]] / expected[["nPerIdent"]]
        expect_identical(object, expected)
    },
    object = list(
        SingleCellExperiment = sce,
        Seurat = seurat
    )
)
