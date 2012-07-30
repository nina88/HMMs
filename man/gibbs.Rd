\name{gibbs}
\alias{gibbs}
\alias{initialise}
\alias{checkpoint}
\alias{initialise_prior}
\alias{initialise_checkpoint}
\alias{log_prior}
\alias{log_likelihood}
\alias{initialise_transition_matrices}
\alias{class_check}
\title{
Gibbs sampling
}
\description{
Gibbs sampling for HMMs
}
\usage{
gibbs(y, iter, prior, r, burnin, thin, checkpoint = NULL)
}
\arguments{
  \item{y}{hmm_fasta object}
  \item{iter}{number of iterations}
  \item{prior}{prior object}
  \item{r}{number of hidden states}
  \item{checkpoint}{checkpoint object if checkpointing is required. Make checkpoint a multiple of thin}
  \item{thin}{what you would like to thin the results by}
  \item{burnin}{burnin size you would like to remove}
}
\details{}
\value{
\item{lambda }{transition matrix for hidden sequence}
\item{P }{array of r transition matrices}
}
\references{}
\author{Nina Wilkinson}
\note{}

\seealso{}
\examples{

}
