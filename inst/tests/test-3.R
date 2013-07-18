test_that("test FB", {
  lambda <- matrix(c(0.3,0.7,0.6,0.4),ncol=2,byrow=TRUE)
  P.1 <- matrix(c(0.3,0.7,0.2,0.8),ncol=2,byrow=TRUE)
  P.2 <- matrix(c(0.5,0.5,0.4,0.6),ncol=2,byrow=TRUE)
  P <- array(0,c(2,2,2))
  P[,,1] <- P.1
  P[,,2] <- P.2
  x=FB(c(1,1,2), lambda, P)
  expect_that(round(x$back,7),  equals(array(c(0.1764706, 0.8235294, 0.4285714, 0.5714286, 0.2379471,  0.7620529, 0.5221843, 0.4778157,0,0,0,0),c(2,2,3))))
  expect_that(round(x$f,7),  equals(matrix(c(0.3, 0.3844221, 0.5683559, 0.7, 0.6155779, 0.4316441),byrow=TRUE,nrow=2)))
  })

test_that("test powerFB", {
  lambda <- matrix(c(0.3,0.7,0.6,0.4),ncol=2,byrow=TRUE)
  P.1 <- matrix(c(0.3,0.7,0.2,0.8),ncol=2,byrow=TRUE)
  P.2 <- matrix(c(0.5,0.5,0.4,0.6),ncol=2,byrow=TRUE)
  P <- array(0,c(2,2,2))
  P[,,1] <- P.1
  P[,,2] <- P.2
  x=FBpower(c(1,1,2), lambda, P, 0.5)
  expect_that(round(x$back,7),  equals(array(c(0.2325672,0.7674328,  0.5147186,  0.4852814, 0.2686629,0.7313371, 0.5625078,  0.4374922, 0,0,0,0),c(2,2,3))))
  expect_that(round(x$f,7),  equals(matrix(c(0.3773705, 0.4235371, 0.5149669, 0.6226295, 0.5764629, 0.4850331),byrow=TRUE,nrow=2)))
})


