## The Seurat wiki describes the changes in v3.0+.
## https://github.com/satijalab/seurat/wiki
suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(basejump) # 0.16.0
    library(Seurat) # 4.1.0
    library(SeuratObject) # 4.0.4
})
load_all(helpers = FALSE)
limit <- structure(1e6, class = "object_size")
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
stopifnot(object.size(object) < limit)
validObject(object)
seurat <- object
use_data(seurat, compress = "xz", overwrite = TRUE)
