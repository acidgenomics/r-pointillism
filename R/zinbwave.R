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
#' Attempt to return stashed zinbwave calcs, or recalculate.
#'
#' zinbwave will calculate `normalizedValues` and `weights` matrices, which will
#' be defined in the [assays()] slot.
#'
#' @family Zero Count Inflation Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams zinbwave::zinbwave
#' @param recalculate `logical`. Force recalculation of weights.
#' @param ... Passthrough arguments to [zinbwave::zinbwave()].
#'
#' @return `SingleCellExperiment`.
#' @export
#'
#' @seealso [zinbwave::zinbwave()].
#'
#' @examples
#' library(bcbioSingleCell)
#' object <- sce_small[seq_len(100L), seq_len(100L)]
#' object <- filterCells(object)
#' zinb <- runZinbwave(object)
#' print(zinb)
#' assayNames(zinb)
runZinbwave <- function(
    object,
    BPPARAM = BiocParallel::SerialParam(),  # nolint
    epsilon = 1e12,
    recalculate = FALSE,
    ...
) {
    stopifnot(is(object, "SingleCellExperiment"))
    # Add an assert check for "Param" class object
    assert_is_a_number(epsilon)
    assert_is_a_bool(recalculate)

    # Early return if weights are already calculated
    if (
        hasZinbwave(object) &&
        !isTRUE(recalculate)
    ) {
        message("Object already has zinbwave weights in assays")
        return(object)
    }

    # Ensure S4 object is coerced to SingleCellExperiment class
    message("Running zinbwave...")
    Y <- as(object, "SingleCellExperiment")  # nolint

    # zinbFit doesn't currently support sparse counts.
    # Ensure they are coerced to a dense matrix.
    # Keep an original copy in case they're sparse, and reslot.
    counts <- counts(Y)
    counts(Y) <- as.matrix(counts(Y))

    # This step can take a long time for large datasets (i.e. hours)
    message(printString(system.time({
        zinb <- zinbwave(
            Y = Y,
            K = 0L,
            BPPARAM = BPPARAM,
            epsilon = epsilon,
            ...
        )
    })))
    stopifnot(hasZinbwave(zinb))

    # Re-slot original sparse counts
    assert_are_identical(dimnames(zinb), dimnames(counts))
    counts(zinb) <- counts

    zinb
}
