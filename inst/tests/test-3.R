test_that("check equilbrium distribution", {
  lambda <- matrix(c(0.999,0.001,0.001,0.999),ncol=2,byrow=TRUE)
  eqdist <- equil(lambda)
  expect_that(eqdist, equals(c(0.5,0.5)))
})

test_that("label-switching", {
  ### when label switching should happen
  lambda <- matrix(c(1,2,3,4),ncol=2,byrow=TRUE)
  P.1 <- matrix(c(1,2,3,4),ncol=2,byrow=TRUE)
  P.2 <- matrix(c(5,6,7,8),ncol=2,byrow=TRUE)
  P <- array(0,c(2,2,2))
  P[,,1] <- P.1
  P[,,2] <- P.2
  x=label_switch(c(1,2,1,2),c(2,1,2,1), 2, lambda, P, 2)
  expect_that(as.vector(x$lambda), equals(c(4,2,3,1)))
  expect_that(x$segment.store2, equals(c(2,1,2,1)))
  expect_that(as.vector(x$P),equals(c(5,7,6,8,1,3,2,4)))
  ### when label switching shouldnt happen
  x=label_switch(c(1,2,1,2),c(1,2,1,2), 2, lambda, P, 2)
  expect_that(x$lambda, equals(lambda))
  expect_that(x$segment.store2, equals(c(1,2,1,2)))
  expect_that(x$P,equals(P))
  #### check happens when having to check max matches - label switching should happen
  x=label_switch(c(1,1,1,1,2,2,2),c(2,2,2,2,1,1,1), 2, lambda, P, 2)
  expect_that(as.vector(x$lambda), equals(c(4,2,3,1)))
  expect_that(x$segment.store2, equals(c(2,2,2,2,1,1,1)))
  expect_that(as.vector(x$P),equals(c(5,7,6,8,1,3,2,4)))
  ### label switching shouldnt happen
  x=label_switch(c(1,1,1,1,2,2,2),c(1,1,1,2,2,2,2), 2, lambda, P, 2)
  expect_that(x$lambda, equals(lambda))
  expect_that(x$segment.store2, equals(c(1,1,1,1,2,2,2)))
  expect_that(x$P,equals(P))
  
})

test_that("check equilbrium distribution", {
  lambda <- matrix(c(0.999,0.001,0.001,0.999),ncol=2,byrow=TRUE)
  eqdist <- equil(lambda)
  expect_that(eqdist, equals(c(0.5,0.5)))
})
