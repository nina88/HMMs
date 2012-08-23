\name{powerpost}
\alias{checkpoint_files_power.R}
\alias{class_check_power.R}
\alias{initialise_power.R}
\alias{initialise_transition_matrices_power.R}
\alias{powerpost}
\title{The power posterior method}
\usage{
  powerpost(N, prior, m, r, y, burnin, checkpoint = NULL)
}
\arguments{
  \item{N}{number of temperature parameters required minus
  1}

  \item{y}{An hmm_fasta object}

  \item{prior}{"prior_parameters" class object}

  \item{m}{number of iterations per temperature parameter}

  \item{r}{number of segment types}

  \item{burnin}{amount of burn in required}

  \item{checkpoint}{"hmm_checkpoint" object or NULL if
  checkpointing not required}
}
\value{
  \item{logpPP }{marginal likelihood for r}
}
\description{
  The power posterior method
}
\keyword{character}

