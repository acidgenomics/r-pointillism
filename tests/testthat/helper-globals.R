data(
    sce,
    package = "acidtest",
    envir = environment()
)
data(
    seurat,
    seurat_all_markers,
    seurat_known_markers,
    package = "pointillism",
    envir = environment()
)

sce <- as(seurat, "SingleCellExperiment")

# nolint start
group_vars <- dplyr::group_vars
# nolint end
