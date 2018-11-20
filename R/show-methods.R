# TODO Improve appearance.
# Remove separator bar and import code from basejump.
# TODO Switch to using `showSlotInfo()` (see bcbioRNASeq.)



#' @name show
#' @author Michael Steinbuagh
#' @inherit methods::show title description details params
NULL



#' @importFrom methods show
#' @aliases NULL
#' @export
methods::show



show.CellCycleMarkers <-  # nolint
    function(object) {
        return(cat("DRAFT METHOD"))

        validObject(object)

        # Include the organism information.
        organism <- metadata(object)[["organism"]]
        release <- metadata(object)[["ensemblRelease"]]

        # Include the gene lengths per phase.
        split <- split(x = object, f = object[["phase"]])
        assert_that(is(split, "SplitDataFrameList"))
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
            paste(class(object), metadata(object)[["version"]]),
            paste0(organism, " (Ensembl ", release, ")"),
            url,
            citation,
            separatorBar,
            paste0(names(genes), " (", lengths, "): ", genes)
        )

        cat(return, sep = "\n")
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("CellCycleMarkers"),
    definition = show.CellCycleMarkers
)



show.CellTypeMarkers <-  # nolint
    function(object) {
        return(cat("DRAFT METHOD"))

        validObject(object)

        # Include the organism information.
        organism <- metadata(object)[["organism"]]
        release <- metadata(object)[["ensemblRelease"]]

        # Include the gene lengths per phase.
        split <- split(x = object, f = object[["cellType"]])
        assert_that(is(split, "SplitDataFrameList"))
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
            paste(class(object), metadata(object)[["version"]]),
            paste0(organism, " (Ensembl ", release, ")"),
            url,
            citation,
            separatorBar,
            paste0(names(genes), " (", lengths, "): ", genes)
        )

        cat(return, sep = "\n")
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("CellTypeMarkers"),
    definition = show.CellTypeMarkers
)



show.KnownMarkers <-  # nolint
    function(object) {
        return(cat("DRAFT METHOD"))
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("KnownMarkers"),
    definition = show.KnownMarkers
)



show.SeuratMarkers <-  # nolint
    function(object) {
        # FIXME
        return(cat("DRAFT METHOD"))

        validObject(object)
        data <- as(object, "DataFrame")

        version <- metadata(object)[["version"]]
        date <- metadata(object)[["date"]]

        gr <- metadata(object)[["GRanges"]]
        organism <- metadata(gr)[["organism"]]
        build <- metadata(gr)[["build"]]
        release <- metadata(gr)[["release"]]

        if ("cluster" %in% colnames(data)) {
            f <- "FindAllMarkers"
            clusters <- levels(data[["cluster"]])
            split <- split(x = data, f = object[["cluster"]])
            assert_that(is(split, "SplitDataFrameList"))
            assert_are_identical(names(split), clusters)
        } else {
            f <- "FindMarkers"
        }

        stop("Not added yet.")
    }



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("SeuratMarkers"),
    definition = show.SeuratMarkers
)



#' @rdname show
#' @export
setMethod(
    f = "show",
    signature = signature("SeuratMarkersPerCluster"),
    definition = show.SeuratMarkers
)
