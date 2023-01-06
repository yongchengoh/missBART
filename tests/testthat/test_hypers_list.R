library(missBART)

test_that("hypers_list produces a list of 4", {
  expect_output(str(tree_list()), "List of 4")
  expect_error(tree_list(prior_alpha = -0.1))
})
