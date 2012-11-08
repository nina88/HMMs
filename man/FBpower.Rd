\name{FBpower}
\alias{checkpoint_files_power}
\alias{FBpower}
\alias{initialise_power}
\title{The forward--backward algorithm for the power posterior method}
\usage{
  FBpower(y, lambda, P, t)
}
\arguments{
  \item{y}{An hmm_fasta object}

  \item{lambda}{hidden sequence transition matrix}

  \item{P}{array of transition matrices for observed
  sequence}

  \item{t}{temperature parameter}
}
\value{
  \item{s }{segmentation} \item{xi }{normalising constant}
  \item{back}{backwards probabilities}
}
\description{
  The forward--backward algorithm for the power posterior
  method
}
\keyword{character}

