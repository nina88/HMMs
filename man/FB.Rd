\name{FB}
\alias{FB}
\title{The forward--backward algorithm}
\usage{
  FB(y, lambda, P)
}
\arguments{
  \item{y}{data from hmm_fasta object}

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

