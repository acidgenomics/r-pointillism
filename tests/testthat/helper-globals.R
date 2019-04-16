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
objects <- list(
    SingleCellExperiment = sce,
    seurat = seurat
)
genes <- head(rownames(seurat))

# nolint start
group_vars <- dplyr::group_vars
# nolint end
