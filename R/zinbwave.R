#' Does Object Contain ZINB-WaVE Weights
#'
#' @family Zero Count Inflation Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `logical`.
#' @export
#'
#' @examples
#' hasZinbwave(sce_small)
hasZinbwave <- function(object) {
    stopifnot(is(object, "SingleCellExperiment"))
    # Require `counts` to always be slotted
    stopifnot("counts" %in% assayNames(object))
    all(c("normalizedValues", "weights") %in% assayNames(object))
}



#' Calculate ZINB-WaVE Weights
#'
#' zinbwave will calculate `normalizedValues` and `weights` matrices, which will
#' be defined in the [assays()] slot. If these values have already been
#' calculated, the object will be returned unmodified, unless
#' `recalculate = TRUE`.
#'
#' @note zinbwave defaults to using [BiocParallel::bpparam()] to register the
#'   number of cores to use (`BPPARAM` argument), but we've found that this
#'   works inconsistently across installations. Currently, we recommend using
#'   [BiocParallel::SerialParam()] by default, which will run in serial, using a
#'   single core. On Linux or macOS, [BiocParallel::MulticoreParam()] should
#'   work to run with multiple cores.
#'
#' @family Zero Count Inflation Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams zinbwave::zinbwave
#' @inheritParams zinbwave::zinbFit
#' @inheritParams zinbwave::zinbModel
#' @param recalculate `logical`. Force recalculation of weights.
#' @param verbose `logical`. Run zinbwave in verbose model (for debugging).
#' @param ... Passthrough arguments to [zinbwave::zinbwave()].
#'
#' @return Modified object (S4 class must extend `SingleCellExperiment`), with
#'   weights added to [assays()].
#' @export
#'
#' @seealso [zinbwave::zinbwave()].
#'
#' @examples
#' Y <- sce_small
#' Y <- bcbioSingleCell::filterCells(Y)
#' # Test using 200 non-zero genes
#' Y <- Y[seq_len(200L), ]
#'
#' # Serial
#' zinb <- runZinbwave(
#'     Y = Y,
#'     BPPARAM = BiocParallel::SerialParam()
#' )
#' print(zinb)
#' assayNames(zinb)
#'
#' # Multicore
#' zinb <- runZinbwave(
#'     Y = Y,
#'     BPPARAM = BiocParallel::MulticoreParam(),
#'     verbose = TRUE
#' )
#' print(zinb)
#' assayNames(zinb)
runZinbwave <- function(
    Y,  # nolint
    K = 0L,  # nolint
    epsilon = 1e12,
    # Use serial by default, for cross-platform compatibility.
    BPPARAM = BiocParallel::SerialParam(),  # nolint
    recalculate = FALSE,
    verbose = FALSE,
    ...
) {
    stopifnot(is(Y, "SingleCellExperiment"))
    # Early return if weights are already calculated.
    if (
        hasZinbwave(Y) &&
        !isTRUE(recalculate)
    ) {
        message("Y already has zinbwave weights in assays")
        return(Y)
    }

    # Assert checks ------------------------------------------------------------
    assert_is_a_number(K)
    assert_is_a_number(epsilon)
    # Require valid BiocParallel bpparam.
    stopifnot(identical(
        attributes(class(BPPARAM))[["package"]],
        "BiocParallel"
    ))
    stopifnot(grepl("Param$", class(BPPARAM)))
    assert_is_a_bool(recalculate)
    assert_is_a_bool(verbose)

    # BiocParallel -------------------------------------------------------------
    # Use a progress bar (only applies to multicore).
    bpprogressbar(BPPARAM) <- TRUE

    # Inform the user whether running in parallel or serial.
    bpparamInfo <- capture.output(BPPARAM)
    message(paste(
        "BiocParallel param registered:",
        paste(bpparamInfo, collapse = "\n"),
        sep = "\n"
    ))

    # Prepare Y object for zinbwave --------------------------------------------
    # We're returnining original object class, with modified assays.
    # Note that `zinbwave()` otherwise returns `SingleCellExperiment`.
    object <- Y

    # Coerce to SingleCellExperiment, for consistency.
    Y <- as(Y, "SingleCellExperiment")

    # zinbFit doesn't currently support sparse counts.
    # Ensure they are coerced to a dense matrix.
    # Keep an original copy in case they're sparse, and reslot before return.
    counts(Y) <- as.matrix(counts(Y))

    # Fit a ZINB regression model ----------------------------------------------
    message("Fitting a ZINB regression model...")
    message(paste(
        "CPU time used:",
        printString(system.time({
            fittedModel <- zinbFit(
                Y = Y,
                K = K,
                epsilon = epsilon,
                BPPARAM = BPPARAM,
                verbose = verbose
            )
        })),
        sep = "\n"
    ))
    message(paste(
        capture.output(print(fittedModel)),
        collapse = "\n"
    ))

    # zinbwave -----------------------------------------------------------------
    message("Running zinbwave...")
    message(paste(
        "CPU time used:",
        printString(system.time({
            zinb <- zinbwave(
                Y = Y,
                fitted_model = fittedModel,
                K = K,
                epsilon = epsilon,
                BPPARAM = BPPARAM,
                verbose = verbose,
                ...
            )
        })),
        sep = "\n"
    ))
    stopifnot(hasZinbwave(zinb))

    # Return -------------------------------------------------------------------
    # Re-slot original raw counts, in case they are sparse.
    counts <- counts(object)
    assert_are_identical(dimnames(zinb), dimnames(counts))
    assays(object) <- assays(zinb)
    counts(object) <- counts
    object
}
