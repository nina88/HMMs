\name{gibbs}
\alias{checkpoint_files.R}
\alias{class_check.R}
\alias{gibbs}
\alias{initialise.R}
\alias{initialise_transition_matrices.R}
\alias{log_likelihood.R}
\alias{log_prior.R}
\title{Gibbs sampling algorithm}
\usage{
  gibbs(y, iter, prior, r, burnin, thin, checkpoint = NULL)
}
\arguments{
  \item{y}{An hmm_fasta object}

  \item{iter}{number of iterations}

  \item{prior}{"prior_parameters" class object}

  \item{r}{number of segment types}

  \item{burnin}{amount of burnin required}

  \item{thin}{amount of thinning wanted}

  \item{checkpoint}{"hmm_checkpoint" object if
  checkpointing wanted or NULL if not wanted}
}
\description{
  Gibbs sampling algorithm
}
\keyword{character}

