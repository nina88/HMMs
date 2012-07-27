\name{gibbs}
\alias{gibbs}
\alias{initialise}
\alias{checkpoint}
\alias{initialise_prior}
\alias{initialise_checkpoint}
\alias{log_prior}
\alias{log_likelihood}
\alias{initialise_transition_matrices}
\title{
Gibbs sampling
}
\description{
Gibbs sampling for HMMs
}
\usage{
gibbs(y, iter, prior, r, checkpoint = NULL)
}
\arguments{
  \item{y}{hmm_fasta object}
  \item{iter}{number of iterations}
  \item{prior}{prior object}
  \item{r}{number of hidden states}
  \item{checkpoint}{checkpoint object if checkpointing is required. If NULL it is carried out every iteration.}
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
