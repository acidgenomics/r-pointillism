#' Show an Object
#'
#' @name show
#' @family S4 Object
#' @author Michael Steinbuagh
#' @importFrom methods show
#' @inherit methods::show
#' @export
#'
#' @examples
#' show(bcb_small)
NULL



.show.CellCycleMarkers <-  # nolint
    function(object) {
        validObject(object)

        # Include the organism information.
        organism <- metadata(object)[["organism"]]
        release <- metadata(object)[["ensemblRelease"]]

        # Include the gene lengths per phase.
        split <- split(x = object, f = object[["phase"]])
        stopifnot(is(split, "SplitDataFrameList"))
        lengths <- nrow(split)
        genes <- sapply(
            X = split,
            FUN = function(x) {
                x <- x[["geneName"]]
                x <- sort(x)
                x <- c(head(x, n = 2L), "...", tail(x, n = 2L))
                paste(x, collapse = " ")
            },
            USE.NAMES = TRUE
        )

        return <- c(
            bold(paste(class(object), metadata(object)[["version"]])),
            bold(paste0(organism, " (Ensembl ", release, ")")),
            url,
            citation,
            separatorBar,
            paste0(bold(names(genes)), " (", lengths, "): ", genes)
        )

        cat(return, sep = "\n")
    }



.show.CellTypeMarkers <-  # nolint
    function(object) {
        validObject(object)

        # Include the organism information.
        organism <- metadata(object)[["organism"]]
        release <- metadata(object)[["ensemblRelease"]]

        # Include the gene lengths per phase.
        split <- split(x = object, f = object[["cellType"]])
        stopifnot(is(split, "SplitDataFrameList"))
        lengths <- nrow(split)
        genes <- sapply(
            X = split,
            FUN = function(x) {
                x <- x[["geneName"]]
                x <- sort(x)
                x <- c(head(x, n = 2L), "...", tail(x, n = 2L))
                paste(x, collapse = " ")
            },
            USE.NAMES = TRUE
        )

        return <- c(
            bold(paste(class(object), metadata(object)[["version"]])),
            bold(paste0(organism, " (Ensembl ", release, ")")),
            url,
            citation,
            separatorBar,
            paste0(bold(names(genes)), " (", lengths, "): ", genes)
        )

        cat(return, sep = "\n")
    }



# FIXME Need to define this.
.show.SeuratMarkers <-  # nolint
    function(object) {
        validObject(object)

        data <- slot(object, name = "data")
        assert_is_all_of(data, "DataFrame")

        # `FindAllMarkers()`
        if ("cluster" %in% colnames(data)) {
            clusters <- levels(data[["cluster"]])
            split <- split(x = data, f = data[["cluster"]])
            stopifnot(is(split, "SplitDataFrameList"))
            assert_are_identical(names(split), clusters)
        }

        stop("Not added yet.")
    }


#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("CellCycleMarkers"),
    definition = .show.CellCycleMarkers
)




#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("CellTypeMarkers"),
    definition = .show.CellTypeMarkers
)
