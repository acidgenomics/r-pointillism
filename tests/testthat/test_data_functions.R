context("Data Functions")



# fetchPCAData =================================================================
test_that("fetchPCAData", {
    x <- fetchPCAData(seurat_small)
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class) %>%
            .[sort(names(.))],
        list(
            centerX = "numeric",
            centerY = "numeric",
            description = "factor",
            ident = "factor",
            index = "factor",
            log10GenesPerUMI = "numeric",
            mitoRatio = "numeric",
            nCoding = "integer",
            nGene = "integer",
            nMito = "integer",
            nUMI = "integer",
            orig.ident = "factor",
            PC1 = "numeric",
            PC2 = "numeric",
            res.0.4 = "character",
            res.0.8 = "character",
            res.1.2 = "character",
            sampleID = "factor",
            sampleName = "factor"
        )
    )
})



# fetchTSNEData ================================================================
test_that("fetchTSNEData", {
    x <- fetchTSNEData(seurat_small)
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class) %>%
            .[sort(names(.))],
        list(
            centerX = "numeric",
            centerY = "numeric",
            description = "factor",
            ident = "factor",
            index = "factor",
            log10GenesPerUMI = "numeric",
            mitoRatio = "numeric",
            nCoding = "integer",
            nGene = "integer",
            nMito = "integer",
            nUMI = "integer",
            orig.ident = "factor",
            res.0.4 = "character",
            res.0.8 = "character",
            res.1.2 = "character",
            sampleID = "factor",
            sampleName = "factor",
            tSNE_1 = "numeric",
            tSNE_2 = "numeric"
        )
    )
})



# fetchTSNEExpressionData ======================================================
test_that("fetchTSNEExpressionData", {
    x <- fetchTSNEExpressionData(
        object = seurat_small,
        genes = head(rownames(seurat_small))
    )
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class) %>%
            .[sort(names(.))],
        list(
            centerX = "numeric",
            centerY = "numeric",
            description = "factor",
            ident = "factor",
            index = "factor",
            log10GenesPerUMI = "numeric",
            mean = "numeric",
            median = "numeric",
            mitoRatio = "numeric",
            nCoding = "integer",
            nGene = "integer",
            nMito = "integer",
            nUMI = "integer",
            orig.ident = "factor",
            res.0.4 = "character",
            res.0.8 = "character",
            res.1.2 = "character",
            sampleID = "factor",
            sampleName = "factor",
            sum = "numeric",
            tSNE_1 = "numeric",
            tSNE_2 = "numeric"
        )
    )
})
