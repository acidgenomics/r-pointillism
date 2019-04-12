data(
    sce,
    package = "acidtest",
    envir = environment()
)
data(
    seurat_small,
    all_markers_small, known_markers_small,
    package = "pointillism",
    envir = environment()
)

# nolint start
group_vars <- dplyr::group_vars
# nolint end

sce <- as(seurat_small, "SingleCellExperiment")
