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
gibbs(y, iter, mu, s, hour, r)
}
\arguments{
  \item{y}{hmm_fasta object}
  \item{iter}{number of iterations}
  \item{mu}{mean}
  \item{s}{standard deviation}
  \item{hour}{checkpoint iteration number}
  \item{r}{number of segments}
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
