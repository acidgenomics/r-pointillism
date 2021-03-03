## Fix for edgeR partial match warning.
options(
    warnPartialMatchAttr = FALSE,
    warnPartialMatchDollar = FALSE
)

data(
    Seurat,
    SingleCellExperiment,
    package = "AcidTest",
    envir = environment()
)
data(
    ## cell_data_set,
    seuratAllMarkers,
    seuratKnownMarkers,
    package = "pointillism",
    envir = environment()
)

## > cds <- cell_data_set
## > rm(cell_data_set)
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
