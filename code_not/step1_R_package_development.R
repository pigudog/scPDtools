# 安装并加载开发包所需包
install.packages("usethis", "devtools", "roxygen2")

library(usethis)
library(devtools)
library(roxygen2)

# 检查
has_devel()
# Your system is ready to build packages!

# create package
create_package("../scPDtools")
# > create_package("../scPDtools")
# ✔ Setting active project to 'D:/program/scPDtools'
# Overwrite pre-existing file 'DESCRIPTION'?
#
#   1: Yup
# 2: Absolutely not
# 3: Not now
#
# Selection: 1
# ✔ Writing 'DESCRIPTION'
# Package: scPDtools
# Title: What the Package Does (One Line, Title Case)
# Version: 0.0.0.9000
# Authors@R (parsed):
#   * First Last <first.last@example.com> [aut, cre] (YOUR-ORCID-ID)
# Description: What the package does (one paragraph).
# License: `use_mit_license()`, `use_gpl3_license()` or friends to
# pick a license
# Encoding: UTF-8
# Roxygen: list(markdown = TRUE)
# RoxygenNote: 7.2.3
# ✔ Writing 'NAMESPACE'
# Overwrite pre-existing file 'scPDtools.Rproj'?
#
#   1: I agree
# 2: Nope
# 3: Not now
#
# Selection: 1
# ✔ Writing 'scPDtools.Rproj'
# ✔ Adding '^scPDtools\\.Rproj$' to '.Rbuildignore'
# ✔ Adding '.Rproj.user' to '.gitignore'
# ✔ Adding '^\\.Rproj\\.user$' to '.Rbuildignore'
# ✔ Opening 'D:/program/scPDtools/' in new RStudio session
# ✔ Setting active project to '<no active project>'
use_data("data/mac_menstrual.rda")
