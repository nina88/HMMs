\name{label_switch}
\alias{initialise_RT}
\alias{label_switch}
\alias{order_of_firsts}
\alias{permut}
\alias{permutations}
\title{Label switching}
\usage{
  label_switch(RT, s, n, lambda, P)
}
\arguments{
  \item{RT}{counts for each simulated sequnce}

  \item{s}{current segmentation sequence}

  \item{n}{length of sequence}

  \item{lambda}{transition matrix for hidden sequence}

  \item{P}{array of transition matrices for observed
  sequence}
}
\value{
  \item{RT}{counts} \item{s}{segmentation after label
  switching corrected for} \item{P }{tranition matrix P
  after label switching corrected for }
  \item{lambda}{trabsition matric lambda after label
  switching corrected for}
}
\description{
  Label switching
}
\keyword{character}

