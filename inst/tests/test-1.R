test_that("trigonometric functions match identies", {
expect_that(sin(pi / 4), equals(1 / sqrt(2)))
expect_that(cos(pi / 4), equals(1 / sqrt(2)))
expect_that(tan(pi / 4), equals(1))
})

test_that("trigonometric functions match identities", {
expect_that(sin(pi / 4), equals(1))
})
