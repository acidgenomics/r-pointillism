## Fix for edgeR partial match warning.
options(
    warnPartialMatchAttr = FALSE,
    warnPartialMatchDollar = FALSE
)

data(
    Seurat,
    SingleCellExperiment,
    cell_data_set,
    package = "acidtest",
    envir = environment()
)
data(
    seuratAllMarkers,
    seuratKnownMarkers,
    package = "pointillism",
    envir = environment()
)

cds <- cell_data_set
rm(cell_data_set)

sce <- SingleCellExperiment
rm(SingleCellExperiment)

seurat <- Seurat
rm(Seurat)

genes <- head(rownames(seurat))

sce <- as(seurat, "SingleCellExperiment")
objects <- list(
    SingleCellExperiment = sce,
    Seurat = seurat
)

## nolint start
group_vars <- dplyr::group_vars
## nolint end
