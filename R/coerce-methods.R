#' Force an object to belong to a class
#'
#' @name coerce
#' @importFrom methods coerce
#' @exportMethod coerce
#' @note Updated 2021-03-03.
#'
#' @seealso
#' - `Seurat::CreateSeuratObject()`.
#' - `Seurat::as.Seurat()`.
#' - `Seurat::as.SingleCellExperiment()`.
#'
#' @examples
#' data(Seurat, SingleCellExperiment, package = "AcidTest")
#'
#' ## SingleCellExperiment to Seurat ====
#' x <- as(SingleCellExperiment, "Seurat")
#' class(x)
#' print(x)
#'
#' ## Seurat to SingleCellExperiment ====
#' x <- as(Seurat, "SingleCellExperiment")
#' print(x)
NULL



## Note that some older saved Seurat objects have row names containing names,
## which are no longer considered valid for a SingleCellExperiment.
##
## For example:
## ENSG00000000003 ENSG00000000005 ENSG00000000419
##        "TSPAN6"          "TNMD"          "DPM1"
## ENSG00000000457 ENSG00000000460 ENSG00000000971
##         "SCYL3"      "C1orf112"           "CFH"
##
## We need to ensure the dimnames themselves are unnamed in the Seurat object
## prior to coercine to SCE, otherwise we will hit a `.validate_names` error.
##
## Inspect current `methods("dimnames")` for Seurat methods:
## > getS3method(f = "dimnames", class = "Seurat")
## > getS3method(f = "dimnames", class = "Assay")
##
## > return(dimnames(x = GetAssay(object = x)))
## > return(dimnames(x = GetAssayData(object = x)))
##
## We need to fix the row names of the primary assay data.
## > data <- GetAssayData(from)
## > class(data)
## ## dgCMatrix
##
## > methods("GetAssayData")
## > getS3method(f = "GetAssayData", class = "Seurat")
## > getS3method(f = "GetAssayData", class = "Assay")
##
## Note that "slot" here defaults to "data" currently.
## > return(slot(object = object, name = slot))
##
## For older objects, `DefaultAssay()` will be "RNA".
## > from@assays$RNA
## [1] "counts"        "data"          "scale.data"
## [4] "key"           "var.features"  "meta.features"
## [7] "misc"
##
## These can contain names and need to be sanitized via `unname()`.
## > rownames(from@assays$RNA@counts)
## > rownames(from@assays$RNA@data)



## Updated 2021-03-03.
`coerce,Seurat,SCE` <-  # nolint
    function(from) {
        validObject(from)
        ## Strip legacy names from row names, if necessary.
        if (hasNames(rownames(GetAssayData(from)))) {
            assay <- DefaultAssay(from)
            assert(
                hasNames(rownames(from@assays[[assay]]@counts)),
                hasNames(rownames(from@assays[[assay]]@data))
            )
            rownames(from@assays[[assay]]@counts) <-
                unname(rownames(from@assays[[assay]]@counts))
            rownames(from@assays[[assay]]@data) <-
                unname(rownames(from@assays[[assay]]@data))
        }
        ## Using the Seurat S3 coercion method here.
        to <- as.SingleCellExperiment(x = from, assay = NULL)
        ## Harden against invalid ident mapping, which can happen when user
        ## reassigns via `Idents<-`.
        ##
        ## Note that Seurat 3.1.3 doesn't currently update "ident" data in
        ## object@meta.data correctly, which is what gets passed to SCE in
        ## default coercion method.
        ##
        ## > Seurat:::Idents.Seurat
        ## > slot(object = object, name = "active.ident")
        ##
        ## > Seurat:::`Idents<-.Seurat`
        ## > slot(object = object, name = "active.ident") <- idents
        ##
        ## Updated on 2020-02-21.
        ## File a bug report with Seurat.
        if (isSubset("ident", colnames(colData(to)))) {
            idents <- Idents(from)
            assert(
                identical(names(idents), rownames(colData(to))),
                is.factor(idents),
                !all(is.na(idents))
            )
            colData(to)[["ident"]] <- unname(idents)
        }
        ## Assays.
        ## > if (isCharacter(assayNames(to))) {
        ## >     assayNames(to) <-
        ## >         camelCase(assayNames(to), strict = TRUE)
        ## > }
        ## Reduced dimensions.
        ## > if (isCharacter(reducedDimNames(to))) {
        ## >     reducedDimNames(to) <-
        ## >         camelCase(reducedDimNames(to), strict = TRUE)
        ## > }
        ## > reducedDims(to) <- lapply(
        ## >     X = reducedDims(to),
        ## >     FUN = function(x) {
        ## >         if (hasColnames(x)) {
        ## >             colnames(x) <-
        ## >                 camelCase(colnames(x), strict = TRUE)
        ## >         }
        ## >         x
        ## >     }
        ## > )
        ## Row and column data.
        rowRanges(to) <- rowRanges(from)
        ## > if (hasColnames(mcols(rowRanges(to)))) {
        ## >     colnames(mcols(rowRanges(to))) <-
        ## >         camelCase(colnames(mcols(rowRanges(to))), strict = TRUE)
        ## > }
        ## > if (hasColnames(colData(to))) {
        ## >     colnames(colData(to)) <-
        ## >         camelCase(colnames(colData(to)), strict = TRUE)
        ## > }
        ## Metadata.
        metadata(to) <- metadata(from)
        metadata(to)[["scaleData"]] <- GetAssayData(from, slot = "scale.data")
        metadata(to)[["variableFeatures"]] <- VariableFeatures(from)
        ## > names(metadata(to)) <-
        ## >     camelCase(names(metadata(to)), strict = TRUE)
        validObject(to)
        to
    }



#' @rdname coerce
#' @name coerce,Seurat,SCE-method
#'
#' @section `Seurat` to `SingleCellExperiment`:
#' S4 coercion support for creating a `SingleCellExperiment` from a `Seurat`
#' class object. The [Seurat FAQ page](https://satijalab.org/Seurat/faq)
#' explains the `Seurat` S4 class structure in detail. Internally, this method
#' improves the basic `Seurat::as.SingleCellExperiment` S3 coercion method,
#' including the `object@scale.data` matrix, and will keep track of stashed
#' `rowRanges` and `metadata` if the `Seurat` object was originally created
#' from a `SingleCellExperiment` (i.e. from the bcbioSingleCell package).
setAs(
    from = "Seurat",
    to = "SingleCellExperiment",
    def = `coerce,Seurat,SCE`
)



## Updated 2019-07-31.
`coerce,Seurat,RSE` <-  # nolint
    function(from) {
        to <- from
        to <- as(to, "SingleCellExperiment")
        to <- as(to, "RangedSummarizedExperiment")
        to
    }



#' @rdname coerce
#' @name coerce,Seurat,RSE-method
#'
#' @section `Seurat` to `RangedSummarizedExperiment`:
#' S4 coercion support for creating a `RangedSummarizedExperiment` from a
#' `Seurat` class object.
setAs(
    from = "Seurat",
    to = "RangedSummarizedExperiment",
    def = `coerce,Seurat,RSE`
)



## Updated 2019-07-31.
`coerce,Seurat,SE` <-  # nolint
    function(from) {
        to <- from
        to <- as(to, "RangedSummarizedExperiment")
        to <- as(to, "SummarizedExperiment")
        to
    }



#' @rdname coerce
#' @name coerce,Seurat,SE-method
#' @section `Seurat` to `SummarizedExperiment`:
#' S4 coercion support for creating a `SummarizedExperiment` from a `Seurat`
#' class object.
setAs(
    from = "Seurat",
    to = "SummarizedExperiment",
    def = `coerce,Seurat,SE`
)



## Updated 2019-07-31.
`coerce,SCE,Seurat` <-  # nolint
    function(from) {
        ## Create the Seurat object. Note that `as.Seurat()` method requires
        ## `logcounts` to be defined in `assays()`, so we're using
        ## `CreateSeuratObject()` here instead.
        to <- CreateSeuratObject(
            counts = counts(from),
            project = "pointillism",
            assay = "RNA",
            min.cells = 0L,
            min.features = 0L,
            names.field = 1L,
            names.delim = "_",
            meta.data = as.data.frame(colData(from))
        )
        ## Check that the dimensions match exactly.
        assert(identical(x = dim(from), y = dim(to)))
        ## Stash rowRanges.
        rowRanges <- rowRanges(from)
        ## Stash metadata.
        metadata <- metadata(from)
        ## Update the session information.
        metadata[["sessionInfo"]] <- session_info()
        ## Seurat v3 still recommends using `misc` slot.
        misc <- list(
            rowRanges = rowRanges,
            metadata = metadata
        )
        misc <- Filter(Negate(is.null), misc)
        slot(to, name = "misc") <- misc
        ## Return.
        to
    }



#' @rdname coerce
#' @name coerce,SCE,Seurat-method
#'
#' @section `SingleCellExperiment` to `Seurat`:
#' Interally `Seurat::CreateSeuratObject` is called without applying any
#' additional filtering cutoffs, since we have already defined them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new `Seurat` class object. Use
#' `convertGenesToSymbols` to convert gene IDs to names (symbols).
setAs(
    from = "SingleCellExperiment",
    to = "Seurat",
    def = `coerce,SCE,Seurat`
)
