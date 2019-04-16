context("plotPCElbow")

# We're testing here to ensure that our seurat matches pbmc_small.

with_parameters_test_that(
    "Seurat", {
        x <- plotPCElbow(object)
        expect_identical(x, seq)
    },
    object = list(
        seurat,
        Seurat::pbmc_small
    ),
    seq = list(
        seq_len(5L),
        seq_len(5L)
    )
)
