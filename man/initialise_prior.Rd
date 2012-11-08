\name{initialise_prior}
\alias{initialise_prior}
\title{Initialising the prior}
\usage{
  initialise_prior(P.mat, mu, s, r, f)
}
\arguments{
  \item{P.mat}{prior value for the observed sequence}

  \item{mu}{prior mean of the hidden sequence}

  \item{s}{prior standard deviation for the hidden
  sequence}

  \item{r}{number of segment types}

  \item{f}{is number of categories in amino acid sequence}
}
\value{
  \item{prior }{"prior_parameters" object}
}
\description{
  Initialising the prior
}
\keyword{character}

