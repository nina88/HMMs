test_that("test convert", {
  expect_that(convert(LETTERS[c(1,3:9,11:14, 16:20,22,23,25)],2),  equals(c(1, 1, 2, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 2, 2, 2, 1, 1, 1)))
  expect_that(convert(LETTERS[c(1,3:9,11:14, 16:20,22,23,25)],3),  equals(c(1, 1, 3, 3, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1)))
  expect_that(convert(LETTERS[c(1,3:9,11:14, 16:20,22,23,25)],4),  equals(c(1, 1, 4, 4, 1, 1, 2, 1, 3, 1, 1, 2, 1, 2, 3, 2, 2, 1, 1, 1)))
  expect_that(convert(LETTERS[c(1,3:9,11:14, 16:20,22,23,25)],5),  throws_error())
  })

test_that("test read_FASTA", {
  x=read_FASTA("inst/extdata/p53.txt",2)
  expect_that(class(x), matches("hmm_fasta"))
})

test_that("test combine_FASTA", {
  x=read_FASTA("inst/extdata/p53.txt",2)
  y=combine_FASTA(x,x)
  expect_that(class(y), matches("hmm_fasta"))
  expect_that(length(y$join), equals(2*length(x$join)-1))
  expect_that(length(y$fasta_seq), equals(2*length(x$fasta_seq)))
})

