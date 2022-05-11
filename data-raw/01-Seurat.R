suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(basejump) # 0.16.0
    library(Seurat) # 4.1.0
    library(SeuratObject) # 4.0.4
})
load_all()
## Restrict object size to 1 MB.
limit <- structure(1e6, class = "object_size")
## Location of example data moved from Seurat to SeuratObject in v4.0.
data(pbmc_small, package = "SeuratObject")
## The Seurat wiki describes the changes in v3.0+.
## https://github.com/satijalab/seurat/wiki
Seurat <- UpdateSeuratObject(pbmc_small)
stopifnot(object.size(Seurat) < limit)
## Add UMAP dimensional reduction to example object.
Seurat <- RunUMAP(
    object = Seurat, dims = seq_len(10L),
    umap.method = "uwot",
    metric = "cosine"
)
## Slot row ranges into the Seurat object.
rowRanges <- makeGRangesFromEnsembl(
    organism = "Homo sapiens",
    level = "genes",
    genomeBuild = "GRCh37",
    ignoreVersion = TRUE
)
rowRanges <- as(rowRanges, "GRanges")
x <- rownames(Seurat)
table <- make.unique(as.character(mcols(rowRanges)[["geneName"]]))
names(rowRanges) <- table
stopifnot(all(x %in% table))
which <- match(x = x, table = table)
rowRanges <- rowRanges[which]
rowRanges <- droplevels2(rowRanges)
stopifnot(object.size(rowRanges) < limit)
rowRanges(Seurat) <- rowRanges
stopifnot(object.size(Seurat) < limit)
validObject(Seurat)
use_data(Seurat, compress = "xz", overwrite = TRUE)
