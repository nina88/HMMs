\name{label_switch}
\alias{label_switch}
\alias{order_of_firsts.R}
\alias{permutations.R}
\alias{permut.R}
\title{Label switching}
\usage{
  label_switch(segment.store2, segment.store1, r, lambda,
    P, f)
}
\arguments{
  \item{segment.store2}{segmentation at iteration i}

  \item{segment.store1}{segmentation at iteration i-1}

  \item{r}{number of segment types}

  \item{lambda}{transition matrix for hidden sequence}

  \item{P}{array of transition matrices for observed
  sequence}

  \item{f}{number of possible states of observed sequence}
}
\value{
  \item{segment.store2 }{segmentation after label switching
  corrected for} \item{P }{tranition matrix P after label
  switching corrected for } \item{lambda}{trabsition matric
  lambda after label switching corrected for}
}
\description{
  Label switching
}
\keyword{character}

