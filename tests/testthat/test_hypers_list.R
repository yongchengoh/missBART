library(missBART)

test_that("hypers_list produces a list of 4", {
  expect_output(str(hypers_list()), "List of 4")
})
