\name{FB}
\alias{backward.R}
\alias{FB}
\alias{forward.R}
\title{The forward--backward algorithm}
\usage{
  FB(y, lambda, P)
}
\arguments{
  \item{y}{An hmm_fasta object}

  \item{lambda}{hidden sequence transition matrix}

  \item{P}{array of transition matrices for observed
  sequence}
}
\value{
  \item{s }{segmentation} \item{xi }{normalising constant}
  \item{back}{backwards probabilities}
}
\description{
  The forward--backward algorithm
}
\keyword{character}

