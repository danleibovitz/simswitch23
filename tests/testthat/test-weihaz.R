
# Set test arguments
len <- 10

# Test negative argument errors ####
test_that("weihaz throws error with a negative scale parameter", {
  expect_error(weihaz(1:len, 1, -1))
})
test_that("weihaz throws error with a negative shape parameter", {
  expect_error(weihaz(1:len, -1, 1))
})
test_that("weihaz throws error with any negative x parameter", {
  expect_error(weihaz(-1:len, 1, 1))
})

# Test argument class errors ####
test_that("weihaz throws error with a character scale parameter", {
  expect_error(weihaz(1:len, 1, "hello"))
})
test_that("weihaz throws error with a character shape parameter", {
  expect_error(weihaz(1:len, "hello", 1))
})
test_that("weihaz throws error with character x parameter", {
  expect_error(weihaz(c("1", "2", "3"), 1, 1))
})

# Test return value always positive ####
test_that("weihaz returns only positive values", {
  expect_gt(min(weihaz(1:len, 1, 1)), 0)
})


# Test the length of the vector ####
test_that("weihaz returns a vector of the correct length", {
  expect_length(weihaz(1:len, 1, 1), length(1:len))
})

