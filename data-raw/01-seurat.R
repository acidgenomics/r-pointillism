## The Seurat wiki describes the changes in v3.0+.
## https://github.com/satijalab/seurat/wiki
## nolint start
suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(basejump)
    library(Seurat)
    library(SeuratObject)
})
## nolint end
# Consider disabling helpers here, if necessary.
load_all()
limit <- structure(1e6L, class = "object_size")
data(pbmc_small, package = "SeuratObject")
object <- pbmc_small
object <- UpdateSeuratObject(object)
stopifnot(object.size(object) < limit)
object <- RunUMAP(object, dims = seq_len(10L))
## Slot row ranges into the Seurat object.
rowRanges <- makeGRangesFromEnsembl(
    organism = "Homo sapiens",
    level = "genes",
    genomeBuild = "GRCh37",
    ignoreVersion = TRUE
)
rowRanges <- as(rowRanges, "GRanges")
rn <- rownames(object)
table <- make.unique(as.character(mcols(rowRanges)[["geneName"]]))
names(rowRanges) <- table
stopifnot(all(rn %in% table))
which <- match(x = rn, table = table)
rowRanges <- rowRanges[which]
rowRanges <- droplevels2(rowRanges)
stopifnot(object.size(rowRanges) < limit)
rowRanges(object) <- rowRanges
stopifnot(
    object.size(object) < limit,
    validObject(object)
)
seurat <- object
use_data(seurat, compress = "xz", overwrite = TRUE)
