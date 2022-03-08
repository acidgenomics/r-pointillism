.formalsList <- list(
    "BPPARAM" = quote(BiocParallel::bpparam()),  # nolint
    "continuousColor" =
        quote(getOption(
            x = "acid.continuous.color",
            default = ggplot2::scale_color_gradient(
                low = "gray75",
                high = "purple"
            )
        )),
    "continuousColorPurpleOrange" =
        quote(getOption(
            x = "acid.continuous.color",
            default = ggplot2::scale_color_gradient2(
                low = "orange",
                mid = "gray75",
                high = "purple",
                midpoint = 0L
            )
        )),
    "discreteColor" =
        quote(getOption(
            x = "acid.discrete.color",
            default = AcidPlots::scale_color_synesthesia_d()
        )),
    "dark" =
        quote(getOption(
            x = "acid.dark",
            default = FALSE
        )),
    "dims" = quote(c(1L, 2L)),
    "direction" = c("both", "up", "down"),
    "expression" = c("mean", "sum"),
    "headerLevel" = 2L,
    "label" =
        quote(getOption(
            x = "acid.label",
            default = TRUE
        )),
    "labelSize" =
        quote(getOption(
            x = "acid.labelSize",
            default = 6L
        )),
    "legend" =
        quote(getOption(
            x = "acid.legend",
            default = TRUE
        )),
    "pointAlpha" =
        quote(getOption(
            x = "acid.pointAlpha",
            default = 0.85
        )),
    "pointSize" =
        quote(getOption(
            x = "acid.pointSize",
            default = 0.75
        )),
    "pointsAsNumbers" =
        quote(getOption(
            x = "acid.pointsAsNumbers",
            default = FALSE
        )),
    ## Alternatively, consider using `1L` here instead.
    "reduction" = "UMAP"
)



## Updated 2021-03-03.
.pkgName <- packageName()
.pkgVersion <- packageVersion(.pkgName)



## Updated 2021-03-03.
.prototypeMetadata <- list(
    "date" = Sys.Date(),
    "packageVersion" = .pkgVersion
)
