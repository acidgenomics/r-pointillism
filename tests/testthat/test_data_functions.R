context("Data Functions")



# fetchReducedDimData ==========================================================
test_that("fetchReducedDimData", {
    x <- fetchReducedDimData(
        object = seurat_small,
        reducedDim = "TSNE"
    )
    expect_is(x, "data.frame")
    expect_identical(
        sort(colnames(x)),
        c(
            "batch",
            "centerX",
            "centerY",
            "expLibSize",
            "group",
            "ident",
            "interestingGroups",
            "log10GenesPerUMI",
            "mitoRatio",
            "nCoding",
            "nGene",
            "nMito",
            "nUMI",
            "orig.ident",
            "res.0.4",
            "res.0.8",
            "res.1.2",
            "sampleID",
            "sampleName",
            "tSNE_1",
            "tSNE_2",
            "x",
            "y"
        )
    )
})



# fetchReducedDimExpressionData ================================================
test_that("fetchReducedDimExpressionData", {
    x <- fetchReducedDimExpressionData(
        object = seurat_small,
        genes = head(rownames(seurat_small)),
        reducedDim = "TSNE"
    )
    expect_is(x, "data.frame")
    expect_identical(
        sort(colnames(x)),
        c(
            "batch",
            "centerX",
            "centerY",
            "expLibSize",
            "group",
            "ident",
            "interestingGroups",
            "log10GenesPerUMI",
            "mean",
            "median",
            "mitoRatio",
            "nCoding",
            "nGene",
            "nMito",
            "nUMI",
            "orig.ident",
            "res.0.4",
            "res.0.8",
            "res.1.2",
            "sampleID",
            "sampleName",
            "sum",
            "tSNE_1",
            "tSNE_2",
            "x",
            "y"
        )
    )
})
