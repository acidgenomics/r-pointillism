# TODO Add MAST, Wilcoxon support?
# TODO Consider adding back zingeR support here.



#' Differential Expression
#'
#' Perform pairwise differential expression across groups of cells by fitting
#' to a zero-inflated negative binomial (ZINB) model using the zinbwave package.
#' Currently supports edgeR and DESeq2 as DE callers.
#'
#' @section zinbwave:
#'
#' This function will run a lot faster if you pre-calculate the ZINB weights
#' using the [runZinbwave()] function, which will stash the weights into the
#' [SingleCellExperiment::weights()] slot of the object. Running zinbwave across
#' the entire set of filtered cells also has greater sensitivity for weight
#' calculations.
#'
#' We are currently using an epsilon setting of `1e12`, as recommended by the
#' ZINB-WaVE integration paper. For more information on the zinbwave package,
#' refer to these materials:
#'
#' - [zinbwave paper](https://doi.org/10.1186/s13059-018-1406-4).
#' - [zinbwave vignette](https://bit.ly/2wtDdpS).
#' - [zinbwave-DESeq2 workflow](https://github.com/mikelove/zinbwave-deseq2).
#'
#' @section edgeR:
#'
#' After estimation of the dispersions and posterior probabilities, the
#' [zinbwave::glmWeightedF()] function is used for statistical inference. This
#' is an adaptation of [edgeR::glmLRT()]. It uses an F-test for which the
#' denominator degrees of freedom are by default adjusted according to the
#' downweighting of excess zeros (`ZI = TRUE`). Also, independent filtering can
#' be performed on the obtained p-values (`independentFiltering = TRUE`). We use
#' the independent filtering strategy that was originally implemented in DESeq2.
#' By default, the average fitted values are used as a filter criterion.
#'
#' @section DESeq2:
#'
#' We're providing preliminary support for DESeq2 as the differential expression
#' caller. It is currently considerably slower for large datasets than edgeR.
#'
#' We're trying to follow the conventions used in DESeq2 for contrasts, defining
#' the name of the factor in the design formula, numerator, and denominator
#' level for the fold change calculations. See [DESeq2::results()] for more
#' information.
#'
#' @section Seurat conventions:
#'
#' Note that Seurat currently uses the convention `cells.1` for the numerator
#' and `cells.2` for the denominator. See [Seurat::DiffExpTest()] for
#' additional information.
#'
#' @note We are currently recommending the ZINB-WaVE method over zingeR, since
#'   it is faster, and has been show to be more sensitive for most single-cell
#'   RNA-seq datasets.
#'
#' @name diffExp
#' @include globals.R
#'
#' @inheritParams general
#' @inheritParams runZinbwave
#' @param numerator `character`. Cells to use in the numerator of the contrast
#'   (e.g. treatment).
#' @param denominator `character`. Cells to use in the denominator of the
#'   contrast (e.g. control).
#' @param caller `string`. Package to use for differential expression calling.
#'   Defaults to `"edgeR"` (faster for large datasets) but `"DESeq2"` is also
#'   supported.
#' @param minCells `scalar integer`. Minimum number of cells required to perform
#'   the differential expression analysis.
#' @param minCellsPerGene `scalar integer`. The minimum number of cells where a
#'   gene is expressed, to pass low expression filtering.
#' @param minCountsPerCell `scalar integer`. Minimum number of counts per cell
#'   for a gene to pass low expression filtering. The number of cells is defined
#'   by `minCellsPerGene`.
#'
#' @return Varies depending on the `caller` argument:
#' - `caller = "edgeR"`: `DEGLRT`.
#' - `caller = "DESeq2"`: Unshrunken `DESeqResults`. Use `lfcShrink()` if
#'   shrunken results are desired.
#'
#' @seealso [Seurat::WhichCells()].
#'
#' @examples
#' data(seurat_small)
#' object <- seurat_small
#'
#' ## Calculate ZINB weights (already stashed).
#' ## seurat_small <- runZinbwave(object)
#'
#' ## Compare expression in cluster 3 relative to 2.
#' ident <- clusterID(object)
#' numerator <- names(ident)[ident == "3"]
#' summary(numerator)
#' denominator <- names(ident)[ident == "2"]
#' summary(denominator)
#'
#' ## edgeR.
#' x <- diffExp(
#'     object = object,
#'     numerator = numerator,
#'     denominator = denominator,
#'     caller = "edgeR"
#' )
#' class(x)
#' summary(x)
#'
#' ## DESeq2.
#' x <- diffExp(
#'     object = object,
#'     numerator = numerator,
#'     denominator = denominator,
#'     caller = "DESeq2"
#' )
#' class(x)
#' summary(x)
NULL



# Internal =====================================================================
.designFormula <- ~group



.underpoweredContrast <- function() {
    warning(paste(
        "Skipping DE.",
        "Underpowered contrast (not enough cells)."
    ), call. = FALSE)
}



# diffExp ======================================================================
.diffExp.SingleCellExperiment <-  # nolint
    function(
        object,
        numerator,
        denominator,
        caller = c("edgeR", "DESeq2"),
        minCells = 2L,  # 10L
        minCellsPerGene = 1L,  # 25L
        minCountsPerCell = 1L,  # 5L
        bpparam  # nolint
    ) {
        # Coerce to standard SCE to ensure fast subsetting.
        object <- as(object, "SingleCellExperiment")

        assert_is_character(numerator)
        assert_is_character(denominator)
        # Early return `NULL` on an imbalanced contrast.
        if (
            length(numerator) < minCells ||
            length(denominator) < minCells
        ) {
            .underpoweredContrast()
            return(NULL)
        }
        assert_are_disjoint_sets(numerator, denominator)
        caller <- match.arg(caller)
        assertIsAnImplicitInteger(minCountsPerCell)
        assertIsAnImplicitInteger(minCellsPerGene)
        assert_all_are_positive(c(minCountsPerCell, minCellsPerGene))
        .assertIsBPPARAM(bpparam)

        assert_is_matrix(.weights(object))
        weightsFun <- metadata(object)[["weights"]]
        assert_is_subset(weightsFun, "zinbwave")

        message(paste(
            "Performing differential expression with",
            paste(weightsFun, caller, sep = "-")
        ))

        # Subset the SCE object to contain the input cells.
        cells <- c(numerator, denominator)
        message(paste(
            paste("Total:", length(cells), "cells"),
            paste("Numerator:", length(numerator), "cells"),
            paste("Denominator:", length(denominator), "cells"),
            sep = "\n"
        ))
        object <- object[, cells]

        # Ensure we're using a sparse matrix to calculate the logical matrix.
        counts <- as(counts(object), "sparseMatrix")

        # Gene filter ----------------------------------------------------------
        message("Applying gene expression low pass filter.")
        message(paste(
            "Requiring at least",
            minCellsPerGene,
            "cells with counts of",
            minCountsPerCell,
            "or more per gene"
        ))

        # Filter the genes based on our expression threshold criteria.
        # Note that this step generates a logical matrix, and will calculate
        # a lot faster when using a sparse matrix (see above).
        genes <- Matrix::rowSums(counts >= minCountsPerCell) >= minCellsPerGene
        genes <- names(genes[genes])
        message(paste(length(genes), "of", nrow(object), "genes passed filter"))

        # Early return NULL if no genes pass.
        if (!length(genes)) {
            warning("No genes passed the low count filter.", call. = FALSE)
            return(NULL)
        }

        # Now subset the object by applying our low pass expression filter.
        object <- object[genes, ]

        # Cell filter ----------------------------------------------------------
        # Inform the user if any cells have been removed.
        trash <- setdiff(cells, colnames(object))
        if (length(trash)) {
            message(paste("Removed", length(trash), "low quality cells"))
        }

        # Resize the numerator and denominator after our QC filters.
        # Early return `NULL` if there are less than n cells in either.
        numerator <- intersect(numerator, colnames(object))
        denominator <- intersect(denominator, colnames(object))
        if (
            length(numerator) < minCells ||
            length(denominator) < minCells
        ) {
            .underpoweredContrast()
            return(NULL)
        }

        # Create a cell factor to define the group.
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
        # Ensure denominator is set as reference.
        group <- relevel(group, ref = "denominator")
        object[["group"]] <- group

        # Set up the design matrix.
        design <- model.matrix(~group)
        metadata(object)[["design"]] <- design

        # Ensure raw counts matrix is dense before running DE.
        counts(object) <- as.matrix(counts(object))

        # Perform differential expression (e.g. `.zinbwave.edgeR`).
        fun <- paste("", weightsFun, caller, sep = ".")
        fun <- get(
            x  = fun,
            envir = asNamespace("pointillism"),
            inherits = FALSE
        )
        assert_is_function(fun)
        fun(object)
    }
formals(.diffExp.SingleCellExperiment)[["bpparam"]] <- bpparam



#' @rdname diffExp
#' @export
setMethod(
    f = "diffExp",
    signature = signature("SingleCellExperiment"),
    definition = .diffExp.SingleCellExperiment
)



#' @rdname diffExp
#' @export
setMethod(
    f = "diffExp",
    signature = signature("seurat"),
    definition = getMethod(
        f = "diffExp",
        signature = signature("SingleCellExperiment")
    )
)



# zinbwave =====================================================================
# Van De Berge and Perraudeau and others have shown the LRT may perform better
# for null hypothesis testing, so we use the LRT. In order to use the Wald test,
# it is recommended to set `useT = TRUE`.
#
# For UMI data, for which the expected counts may be very low, the likelihood
# ratio test implemented in nbinomLRT should be used.
#
# DESeq2 supports `weights` in assays automatically.
.zinbwave.DESeq2 <- function(object) {  # nolint
    # TODO Switch to using `design()` generic.
    .assertHasDesignFormula(object)
    # DESeq2 -------------------------------------------------------------------
    message("Running DESeq2.")
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
        # We already performed low count filtering.
        res <- results(dds, independentFiltering = FALSE)
    })))
    res
}



.zinbwave.edgeR <- function(object) {  # nolint
    # TODO Switch to using `design()` generic.
    .assertHasDesignFormula(object)
    # edgeR --------------------------------------------------------------------
    message("Running edgeR.")
    # Coerce to dense matrix.
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
        # We already performed low count filtering.
        lrt <- glmWeightedF(
            glmfit = fit,
            coef = 2L,
            independentFiltering = FALSE
        )
    })))
    lrt
}
