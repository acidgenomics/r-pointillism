rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "warning")
BiocCheck::BiocCheck(`quit-with-status` = TRUE)
lintr::lint_package()
