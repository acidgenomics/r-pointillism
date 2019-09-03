globalVariables(".")

url <- "https://steinbaugh.com/pointillism"
citation <- "citation(\"pointillism\")"

## DR marker default color palettes.
darkMarkerColors <-
    ggplot2::scale_colour_viridis_c(option = "plasma")
lightMarkerColors <-
    ggplot2::scale_colour_gradient(low = "gray90", high = "red")

BPPARAM <- quote(BiocParallel::bpparam())  # nolint

continuousColor <-
    quote(getOption(
        "acid.continuous.color",
        default = ggplot2::scale_colour_gradient(
            low = "gray75",
            high = "purple"
        )
    ))
continuousColorPurpleOrange <-
    quote(getOption(
        "acid.continuous.color",
        default = ggplot2::scale_colour_gradient2(
            low = "orange",
            mid = "gray75",
            high = "purple",
            midpoint = 0L
        )
    ))
discreteColor <-
    quote(getOption("acid.discrete.color", default = NULL))

## FIXME Work on improving consistency here with acidplots.

dark <- quote(getOption("acid.dark", default = FALSE))
dims <- quote(c(1L, 2L))
direction <- c("up", "down", "both")
expression <- c("mean", "sum")
headerLevel <- 2L
label <- quote(getOption("acid.label", default = TRUE))
labelSize <- quote(getOption("acid.labelSize", default = 6L))
legend <- quote(getOption("acid.legend", default = TRUE))
pointAlpha <- quote(getOption("acid.pointAlpha", default = 0.85))
pointSize <- quote(getOption("acid.pointSize", default = 0.75))
pointsAsNumbers <-
    quote(getOption("acid.pointsAsNumbers", default = FALSE))
reduction <- "UMAP"  # 1L

## Set default ggplot2 theme.
if (isTRUE(getOption("acid.dark"))) {
    theme <- acidplots::acid_theme_dark
} else {
    theme <- acidplots::acid_theme_light
}
theme <- theme(base_size = 14L, legend_position = "right")
ggplot2::theme_set(new = theme)
