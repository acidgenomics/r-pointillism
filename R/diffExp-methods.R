#' Differential Expression
#'
#' Perform pairwise differential expression across groups of cells by fitting
#' to a zero-inflated negative binomial (ZINB) model using the zinbwave package.
#' Currently supports edgeR and DESeq2 as DE callers.
#'
#' This step will run a lot faster if you pre-calculate the ZINB weights using
#' the [runZinbwave()] function, which will stash the weights into the
#' [assays()] slot of the object. Running zinbwave across the entire set of
#' filtered cells also has greater sensitivity for weight calculations.
#'
#' @section zinbwave:
#' We are currently using an epsilon setting of `1e12`, as recommended by the
#' ZINB-WaVE integration paper. For more information on the zinbwave package,
#' refer to these materials:
#'
#' - [zinbwave paper](https://doi.org/10.1186/s13059-018-1406-4).
#' - [zinbwave vignette](https://bit.ly/2wtDdpS).
#' - [zinbwave-DESeq2 workflow](https://github.com/mikelove/zinbwave-deseq2).
#'
#' @section edgeR:
#' After estimation of the dispersions and posterior probabilities, the
#' `glmWeightedF()` function is used for statistical inference. This is an
#' adapted function from the `glmLRT()` function of edgeR. It uses an F-test for
#' which the denominator degrees of freedom are by default adjusted according to
#' the downweighting of excess zeros (`ZI = TRUE`). Also, independent filtering
#' can be performed on the obtained p-values (`independentFiltering = TRUE`). We
#' use the independent filtering strategy that was originally implemented in
#' DESeq2. By default, the average fitted values are used as a filter criterion.
#'
#' @section DESeq2:
#' We're providing preliminary support for DESeq2 as the differential expression
#' caller. It is currently considerably slower for large datasets than edgeR.
#'
#' We're trying to follow the conventions used in DESeq2 for contrasts, defining
#' the name of the factor in the design formula, numerator, and denominator
#' level for the fold change calculations. See [DESeq2::results()] for more
#' information.
#'
#' @section Seurat conventions:
#' Note that Seurat currently uses the convention `cells.1` for the numerator
#' and `cells.2` for the denominator. See [Seurat::DiffExpTest()] for
#' additional information.
#'
#' @note We are currently recommending the ZINB-WaVE method over zingeR, since
#'   it is faster, and has been show to be more sensitive for most single-cell
#'   RNA-seq datasets.
#'
#' @name diffExp
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param numerator `character`. Cells to use in the numerator of the contrast
#'   (e.g. treatment).
#' @param denominator `character`. Cells to use in the denominator of the
#'   contrast (e.g. control).
#' @param caller `string`. Package to use for differential expression calling.
#'   Defaults to `"edgeR"` (faster for large datasets) but `"DESeq2"` is also
#'   supported.
#' @param minCellsPerGene `scalar integer`. The minimum number of cells where a
#'   gene is expressed, to pass low expression filtering. Set to `0` to disable
#'   (*not recommended*).
#' @param minCountsPerCell `scalar integer`. Minimum number of counts per cell
#'   for a gene to pass low expression filtering. The number of cells is defined
#'   by `minCellsPerGene`. Set to `0` to disable (*not recommended*).
#'
#' @return Varies depending on the `caller` argument:
#' - `caller = "edgeR"`: `DEGLRT`.
#' - `caller = "DESeq2"`: Unshrunken `DESeqResults`. Use `lfcShrink()` if
#'   shrunken results are desired.
#'
#' @seealso [Seurat::WhichCells()].
#'
#' @examples
#' library(bcbioSingleCell)
#' object <- sce_small[seq_len(100L), seq_len(100L)]
#' object <- filterCells(object)
#'
#' numerator <- colnames(object)[object$group == "group2"]
#' glimpse(numerator)
#'
#' denominator <- colnames(object)[object$group == "group1"]
#' glimpse(denominator)
#'
#' # edgeR
#' x <- diffExp(
#'     object = object,
#'     numerator = numerator,
#'     denominator = denominator,
#'     caller = "edgeR"
#' )
#' class(x)
#' summary(x)
#'
#' # DESeq2
#' x <- diffExp(
#'     object = object,
#'     numerator = numerator,
#'     denominator = denominator,
#'     caller = "DESeq2"
#' )
#' class(x)
#' summary(x)
NULL



.designFormula <- ~group



# Van De Berge and Perraudeau and others have shown the LRT may perform better
# for null hypothesis testing, so we use the LRT. In order to use the Wald test,
# it is recommended to set `useT = TRUE`.
#
# For UMI data, for which the expected counts may be very low, the likelihood
# ratio test implemented in nbinomLRT should be used.
#
# DESeq2 supports `weights` in assays automatically.
.zinbwave.DESeq2 <- function(object) {  # nolint
    .assertHasZinbwave(object)
    .assertHasDesignFormula(object)
    # DESeq2 -------------------------------------------------------------------
    message("Running DESeq2...")
    message(printString(system.time({
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
            minReplicatesForReplace = Inf
        )
        # We already performed low count filtering
        res <- results(dds, independentFiltering = FALSE)
    })))
    res
}



.zinbwave.edgeR <- function(object) {  # nolint
    .assertHasZinbwave(object)
    .assertHasDesignFormula(object)
    # edgeR --------------------------------------------------------------------
    message("Running edgeR...")
    counts <- as.matrix(counts(object))
    weights <- assay(object, "weights")
    assert_is_matrix(weights)
    design <- metadata(object)[["design"]]
    assert_is_matrix(design)
    group <- object[["group"]]
    assert_is_factor(group)
    message(printString(system.time({
        dge <- DGEList(counts, group = group)
        dge <- calcNormFactors(dge)
        dge[["weights"]] <- weights
        dge <- estimateDisp(dge, design = design)
        fit <- glmFit(dge, design = design)
        # We already performed low count filtering
        lrt <- glmWeightedF(
            glmfit = fit,
            coef = 2L,
            independentFiltering = FALSE
        )
    })))
    lrt
}



#' @rdname diffExp
#' @export
setMethod(
    "diffExp",
    signature("SingleCellExperiment"),
    function(
        object,
        numerator,
        denominator,
        caller = c("edgeR", "DESeq2"),
        minCellsPerGene = 25L,
        minCountsPerCell = 5L
    ) {
        object <- as(object, "SingleCellExperiment")
        assert_is_character(numerator)
        assert_is_character(denominator)
        assert_are_disjoint_sets(numerator, denominator)
        caller <- match.arg(caller)
        # Consider adding zingeR support back once it's on Bioconductor
        zeroWeights <- "zinbwave"
        assertIsAnImplicitInteger(minCountsPerCell)
        assertIsAnImplicitInteger(minCellsPerGene)
        assert_all_are_positive(c(minCountsPerCell, minCellsPerGene))

        message(paste(
            "Performing differential expression with",
            paste(zeroWeights, caller, sep = "-")
        ))

        # Subset the SCE object to contain the desired cells
        cells <- c(numerator, denominator)
        object <- object[, cells]

        message(paste(
            "Applying low gene count filter",
            paste(
                "Requiring at least",
                minCellsPerGene,
                "cells with counts of",
                minCountsPerCell,
                "or more per gene"
            ),
            sep = "\n"
        ))
        genes <- rowSums(counts(object) >= minCountsPerCell) >= minCellsPerGene
        # Early return NULL if no genes pass
        if (!any(genes)) {
            warning("No genes passed the low count filter")
            return(NULL)
        }
        object <- object[genes, ]

        # Create a cell factor to define the group for `DGEList()`
        numeratorFactor <- replicate(
            n = length(numerator),
            expr = "numerator"
        ) %>%
            as.factor() %>%
            set_names(numerator)
        denominatorFactor <- replicate(
            n = length(denominator),
            expr = "denominator"
        ) %>%
            as.factor() %>%
            set_names(denominator)
        group <- factor(c(
            as.character(numeratorFactor),
            as.character(denominatorFactor)
        ))
        names(group) <- c(names(numeratorFactor), names(denominatorFactor))
        # Ensure denominator is set as reference
        group <- relevel(group, ref = "denominator")
        object[["group"]] <- group

        # Set up the design matrix
        design <- model.matrix(~group)
        metadata(object)[["design"]] <- design

        # Calculate the weights (e.g. `runZinbwave`)
        weightsFunction <- get(paste0("run", upperCamel(zeroWeights)))
        object <- weightsFunction(Y = object, BPPARAM = SerialParam())

        # Ensure raw counts matrix is dense
        counts(object) <- as.matrix(counts(object))

        # Perform differential expression (e.g. `.zinbwave.edgeR`)
        callerFunction <- get(paste("", zeroWeights, caller, sep = "."))
        object <- callerFunction(object)

        object
    }
)



#' @rdname diffExp
#' @export
setMethod(
    "diffExp",
    signature("seurat"),
    getMethod("diffExp", "SingleCellExperiment")
)
