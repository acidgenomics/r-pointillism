.stashAttributes <- function(object, release) {
    assert_is_a_number(release)
    # Stash the pointillism version.
    attr(object, "version") <- packageVersion("pointillism")
    # Ensure we're stashing the Ensembl release.
    attr(object, "ensemblRelease") <- release
    # Stash the date.
    attr(object, "date") <- Sys.Date()
    object
}
