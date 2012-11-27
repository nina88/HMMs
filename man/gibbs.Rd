\name{gibbs}
\alias{checkpoint_files}
\alias{class_check}
\alias{gibbs}
\alias{initialise}
\alias{initialise_transition_matrices}
\alias{log_likelihood}
\alias{log_prior}
\title{Gibbs sampling algorithm}
\usage{
  gibbs(y, iter, prior, burnin, thin, checkpoint = NULL)
}
\arguments{
  \item{y}{An hmm_fasta object}

  \item{iter}{number of iterations}

  \item{prior}{"prior_parameters" class object}

  \item{burnin}{amount of burnin required}

  \item{thin}{amount of thinning wanted}

  \item{checkpoint}{"hmm_checkpoint" object if
  checkpointing wanted or NULL if not wanted}
}
\description{
  Gibbs sampling algorithm
}
\keyword{character}

