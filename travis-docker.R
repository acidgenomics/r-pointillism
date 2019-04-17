# Avoid time zone error due to `timedatectl` being installed on image.
#
# Warning message:
# In system("timedatectl", intern = TRUE) :
#   running command 'timedatectl' had status 1
#
# https://github.com/rocker-org/rocker-versioned/issues/89
Sys.setenv(TZ = "America/New_York")

utils::sessionInfo()
sessioninfo::session_info()

rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "error")
BiocCheck::BiocCheck(`quit-with-status` = FALSE)

lintr::lint_package()
