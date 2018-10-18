# Require zero weights to be calculated, showing a more informative error.
.weights <- function(object) {
    tryCatch(
        expr = weights(object),
        error = function(e) {
            stop(paste(
                "Object does not contain zero weights in `weights()` slot.",
                "`runZinbwave()` is recommended by default.",
                sep = "\n"
            ), call. = FALSE)
        }
    )
}
