library(testthat)
library(CellMembrane)

test_check("CellMembrane", reporter = "progress")

# if (requireNamespace("lintr", quietly = TRUE)) {
#     context("lints")
#     test_that("Package Style", {
#         lintr::expect_lint_free()
#     })
# }