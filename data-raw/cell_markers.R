# Cell Cycle and Cell Type Markers
# 2018-09-19
# This code is derived from:
# - Tirosh et al, 2015
# - http://satijalab.org/seurat/cell_cycle_vignette.html

library(devtools)
library(googlesheets)
library(tidyverse)

# Must be interactive, requiring Google Sheets authentication.
stopifnot(interactive())

# Ensembl release version.
release <- 92L

# Here we're matching the stored Ensembl identifiers (`geneID`) using
# ensembldb to obtain the latest symbol names from Ensembl.

# Generate an OAuth token on an R server using httr.
# https://support.rstudio.com/hc/en-us/articles/217952868-Generating-OAuth-tokens-for-a-server-using-httr
# Otherwise you will run into a localhost port 1410 error.
# Easiest way to fix this is to set `options(httr_oob_default = TRUE)`.
# You'll get a code that you have to paste back into an R prompt.

# Allow tidyverse to access Google Sheets.
gs_ls()

# Cell cycle markers ===========================================================
# Download the Google sheet (gs).
gs <- gs_key("1qA5ktYeimNGpZF1UPSQZATbpzEqgyxN6daoMOjv6YYw")

# Get a list of the worksheets (ws).
ws <- gs_ws_ls(gs)
print(ws)

cell_cycle_markers <- lapply(ws, function(ws) {
    data <- gs %>%
        gs_read(ws = ws) %>%
        select(phase, geneID) %>%
        mutate(
            geneName = convertGenesToSymbols(
                geneID,
                organism = ws,
                release = release
            )
        ) %>%
        group_by(phase) %>%
        arrange(geneID, .by_group = TRUE)
    data <- .stashAttributes(data, release = release)
    new("CellCycleMarkers", data)
})
names(cell_cycle_markers) <- camel(ws)

# Cell type markers ============================================================
# Download the Google sheet (gs).
gs <- gs_key("1vGNU2CCxpaoTCLvzOxK1hf5gjULrf2-CpgCp9bOfGJ0")

# Get a list of the worksheets (ws).
ws <- gs_ws_ls(gs) %>%
    # Remove internal worksheets prefixed with "_".
    .[!str_detect(., "^_")]
print(ws)

cell_type_markers <- lapply(ws, function(ws) {
    data <- gs %>%
        gs_read(ws = ws) %>%
        select(cellType, geneID) %>%
        mutate(
            geneName = convertGenesToSymbols(
                geneID,
                organism = ws,
                release = release
            )
        ) %>%
        group_by(cellType) %>%
        arrange(geneID, .by_group = TRUE)
    data <- .stashAttributes(data, release = release)
    new("CellTypeMarkers", data)
})
names(cell_type_markers) <- camel(ws)

# Save R data ==================================================================
devtools::use_data(
    cell_cycle_markers,
    cell_type_markers,
    compress = "xz",
    overwrite = TRUE
)
