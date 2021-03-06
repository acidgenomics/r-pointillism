#' @name diffExp
#' @inherit AcidGenerics::diffExp
#'
#' @note We are no longer recommending the use of software that attempts to
#'   mitigate zero count inflation (e.g. zinbwave, zingeR) for UMI droplet-based
#'   single cell RNA-seq data. Simply model the counts directly.
#' @note Updated 2020-02-21.
#'
#' @details
#' Perform pairwise differential expression across groups of cells. Currently
#' supports edgeR and DESeq2 as DE callers.
#'
#' @section DESeq2:
#'
#' We're providing preliminary support for DESeq2 as the differential expression
#' caller. It is currently considerably slower for large datasets than edgeR.
#'
#' We're trying to follow the conventions used in DESeq2 for contrasts, defining
#' the name of the factor in the design formula, numerator, and denominator
#' level for the fold change calculations. See [DESeq2::results()] for details.
#'
#' Van de Berge and Perraudeau and others have shown the LRT may perform better
#' for null hypothesis testing, so we use the LRT. In order to use the Wald
#' test, it is recommended to set `useT = TRUE` (*not currently in use*).
#'
#' For UMI data, for which the expected counts may be very low, the likelihood
#' ratio test implemented in `nbinomLRT()` should be used.
#'
#' Note that DESeq2 supports `weights()` values automatically, if slotted using
#' zinbwave (which is no longer recommended for droplet scRNA-seq).
#'
#' @section edgeR:
#'
#' The LRT has been shown to perform better for null hypothesis testing with
#' droplet scRNA-seq data. Here we are using [edgeR::glmLRT()] internally.
#'
#' edgeR is currently significantly faster than DESeq2 for large datasets.
#'
#' @section Seurat conventions:
#'
#' Note that Seurat currently uses the convention `cells.1` for the numerator
#' and `cells.2` for the denominator. See [Seurat::FindMarkers()] for details.
#'
#' @inheritParams AcidRoxygen::params
#' @param numerator `character`.
#'   Cells to use in the numerator of the contrast (e.g. treatment).
#' @param denominator `character`.
#'   Cells to use in the denominator of the contrast (e.g. control).
#' @param caller `character(1)`.
#'   Package to use for differential expression calling. Defaults to `"edgeR"`
#'   (faster for large datasets) but `"DESeq2"` is also supported.
#' @param minCells `integer(1)`.
#'   Minimum number of cells required to perform the differential expression
#'   analysis.
#' @param minCellsPerGene `integer(1)`.
#'   The minimum number of cells where a gene is expressed, to pass low
#'   expression filtering.
#' @param minCountsPerCell `integer(1)`.
#'   Minimum number of counts per cell for a gene to pass low expression
#'   filtering. The number of cells is defined by `minCellsPerGene`.
#' @param BPPARAM `bpparamClass`.
#'   Back-end method to be used for computations.
#'   See [`bpparam`][BiocParallel::bpparam] for details.
#'   Currently only used by DESeq2 but not edgeR for calculations here.
#' @param ... Additional arguments.
#'
#' @return Varies depending on the `caller` argument:
#'
#' - `caller = "edgeR"`: `DEGLRT`.
#' - `caller = "DESeq2"`: Unshrunken `DESeqResults`.
#'   Apply [DESeq2::lfcShrink()] if shrunken results are desired.
#'
#' @seealso [Seurat::WhichCells()].
#'
#' @examples
#' data(Seurat, package = "AcidTest")
#' object <- Seurat
#'
#' ## Compare expression in cluster 2 relative to 1.
#' clusters <- clusters(object)
#' numerator <- names(clusters)[clusters == "2"]
#' summary(numerator)
#' denominator <- names(clusters)[clusters == "1"]
#' summary(denominator)
#'
#' ## edgeR ====
#' ## > x <- diffExp(
#' ## >     object = object,
#' ## >     numerator = numerator,
#' ## >     denominator = denominator,
#' ## >     caller = "edgeR"
#' ## > )
#' ## > class(x)
#' ## > summary(x)
#'
#' ## DESeq2 ====
#' ## This will warn about weights with the minimal example.
#' ## > x <- diffExp(
#' ## >     object = object,
#' ## >     numerator = numerator,
#' ## >     denominator = denominator,
#' ## >     caller = "DESeq2"
#' ## > )
#' ## > class(x)
#' ## > summary(x)
NULL



.designFormula <- ~group



## DESeq2 is slow for large datasets.
##
## - `reduced`: For `test = "LRT"`, a reduced formula to compare against.
## - `sfType`: Use "poscounts" instead of "ratio" here because we're
##   expecting genes with zero counts.
##   See `DESeq2::estimateSizeFactors()` for details.
## - `minmu`: Set a lower threshold than the default 0.5, as recommended
##   in Mike Love's zinbwave-DESeq2 vignette.
##
## Updated 2020-01-30.
.diffExp.DESeq2 <- function(object, BPPARAM) {  # nolint
    assert(.hasDesignFormula(object))
    alert("Running {.pkg DESeq2}.")
    dds <- DESeqDataSet(
        se = object,
        design = .designFormula
    )
    dds <- DESeq(
        object = dds,
        test = "LRT",
        reduced = ~ 1L,
        sfType = "poscounts",
        minmu = 1e-6,
        minReplicatesForReplace = Inf,
        BPPARAM = BPPARAM
    )
    ## We have already performed low count filtering.
    res <- results(
        object = dds,
        independentFiltering = FALSE,
        BPPARAM = BPPARAM
    )
    res
}



## edgeR is much faster than DESeq2 for large datasets.
##
## Note that zinbwave recommends `glmWeightedF()`, which recycles an old version
## of the `glmLRT()` method, that allows an F-test with adjusted denominator
## degrees of freedom, to account for the downweighting in the zero-inflation
## model (which no longer applies here).
##
## Updated 2020-01-30.
.diffExp.edgeR <- function(object) {  # nolint
    assert(.hasDesignFormula(object))
    alert("Running {.pkg edgeR}.")
    ## Ensure sparseMatrix gets coerced to dense matrix.
    counts <- as.matrix(counts(object))
    design <- metadata(object)[["design"]]
    assert(is.matrix(design))
    group <- object[["group"]]
    assert(is.factor(group))
    dge <- DGEList(counts, group = group)
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design = design)
    fit <- glmFit(dge, design = design)
    lrt <- glmLRT(glmfit = fit, coef = 2L)
    lrt
}



.underpoweredContrast <- function() {
    warning("Skipping DE. Underpowered contrast (not enough cells).")
}



## Updated 2021-03-02.
`diffExp,SingleCellExperiment` <-  # nolint
    function(
        object,
        numerator,
        denominator,
        caller = c("edgeR", "DESeq2"),
        minCells = 2L,  # 10L
        minCellsPerGene = 1L,  # 25L
        minCountsPerCell = 1L,  # 5L
        BPPARAM  # nolint
    ) {
        ## Coerce to standard SCE to ensure fast subsetting.
        object <- as(object, "SingleCellExperiment")
        assert(
            is.character(numerator),
            is.character(denominator)
        )
        ## Early return `NULL` on an imbalanced contrast.
        if (
            length(numerator) < minCells ||
            length(denominator) < minCells
        ) {
            .underpoweredContrast()
            return(NULL)
        }
        assert(
            areDisjointSets(numerator, denominator),
            isInt(minCountsPerCell),
            isInt(minCellsPerGene),
            allArePositive(c(minCountsPerCell, minCellsPerGene)),
            .isBPPARAM(BPPARAM)
        )
        caller <- match.arg(caller)
        alert(sprintf(
            "Performing differential expression with {.pkg %s}.",
            caller
        ))
        ## Subset the SCE object to contain the input cells.
        cells <- c(numerator, denominator)
        ul(c(
            sprintf(
                "Total: %d %s",
                length(cells),
                ngettext(length(cells), "cell", "cells")
            ),
            sprintf(
                "Numerator: %d %s",
                length(numerator),
                ngettext(length(numerator), "cell", "cells")
            ),
            sprintf(
                "Denominator: %d %s",
                length(denominator),
                ngettext(length(denominator), "cell", "cells")
            )
        ))
        object <- object[, cells, drop = FALSE]
        ## Ensure we're using a sparse matrix to calculate the logical matrix.
        counts <- counts(object)
        counts <- as(counts, "sparseMatrix")
        ## Gene filter ---------------------------------------------------------
        alert("Applying gene expression low pass filter.")
        alertInfo(sprintf(
            "Requiring at least %d %s with counts of %d or more per gene.",
            minCellsPerGene,
            ngettext(
                n = minCellsPerGene,
                msg1 = "cell",
                msg2 = "cells"
            ),
            minCountsPerCell
        ))
        ## Filter the genes based on our expression threshold criteria.
        ## Note that this step generates a logical matrix, and will calculate
        ## a lot faster when using a sparse matrix (see above).
        genes <- rowSums(counts >= minCountsPerCell) >= minCellsPerGene
        genes <- names(genes[genes])
        alertInfo(sprintf(
            "%d of %d %s passed filter.",
            length(genes),
            nrow(object),
            ngettext(
                n = nrow(object),
                msg1 = "gene",
                msg2 = "genes"
            )
        ))
        ## Early return NULL if no genes pass.
        if (!hasLength(genes)) {
            alertWarning("No genes passed the low count filter.")
            return(NULL)
        }
        ## Now subset the object by applying our low pass expression filter.
        object <- object[genes, , drop = FALSE]
        ## Cell filter ---------------------------------------------------------
        alert("Applying cell low pass filter.")
        ## Inform the user if any cells have been removed.
        trash <- setdiff(cells, colnames(object))
        if (hasLength(trash)) {
            alertWarning(sprintf(
                "Removed %d low quality %s.",
                length(trash),
                ngettext(
                    n = length(trash),
                    msg1 = "cell",
                    msg2 = "cells"
                )
            ))
        }
        ## Resize the numerator and denominator after our QC filters.
        ## Early return `NULL` if there are less than n cells in either.
        numerator <- intersect(numerator, colnames(object))
        denominator <- intersect(denominator, colnames(object))
        if (
            length(numerator) < minCells ||
            length(denominator) < minCells
        ) {
            .underpoweredContrast()
            return(NULL)
        }
        ## Design formula ------------------------------------------------------
        ## Create a cell factor to define the group.
        numeratorFactor <- as.factor(replicate(
            n = length(numerator),
            expr = "numerator"
        ))
        names(numeratorFactor) <- numerator
        denominatorFactor <- as.factor(replicate(
            n = length(denominator),
            expr = "denominator"
        ))
        names(denominatorFactor) <- denominator
        group <- factor(c(
            as.character(numeratorFactor),
            as.character(denominatorFactor)
        ))
        names(group) <- c(names(numeratorFactor), names(denominatorFactor))
        ## Ensure denominator is set as reference.
        group <- relevel(group, ref = "denominator")
        object[["group"]] <- group
        ## Set up the design matrix.
        design <- model.matrix(~group)
        metadata(object)[["design"]] <- design
        ## Run differential expression -----------------------------------------
        ## Ensure count matrix is dense before running DE.
        counts(object) <- as.matrix(counts(object))
        ## Perform differential expression.
        fun <- get(
            x  = paste0(".diffExp.", caller),
            envir = asNamespace("pointillism"),
            inherits = FALSE
        )
        assert(is.function(fun))
        fun(object)
    }

args <- "BPPARAM"
formals(`diffExp,SingleCellExperiment`)[args] <-
    .formalsList[args]
rm(args)



#' @rdname diffExp
#' @export
setMethod(
    f = "diffExp",
    signature = signature("SingleCellExperiment"),
    definition = `diffExp,SingleCellExperiment`
)



## Updated 2020-02-21.
`diffExp,Seurat` <-  # nolint
    function(object, ...) {
        diffExp(object = as(object, "SingleCellExperiment"), ...)
    }



#' @rdname diffExp
#' @export
setMethod(
    f = "diffExp",
    signature = signature("Seurat"),
    definition = `diffExp,Seurat`
)
