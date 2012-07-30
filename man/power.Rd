\name{powerpost}
\alias{powerpost}
\alias{FBpower}
\alias{checkpoint_files_power}
\alias{initialise_power}
\alias{initialise_checkpoint_power}
\alias{initialise_transition_matrices_power}
\alias{class_check_power}
\title{
Power Posterior method
}
\description{
Power posterior method for finding the marginal likelihood
}
\usage{
powerpost(N, prior, m, r, y, burnin, checkpoint=NULL)
}
\arguments{
  \item{N}{number of values of t minus 1}
  \item{m}{number of iterations per value of t}
  \item{prior}{prior object}
  \item{r}{number of hidden states}
  \item{checkpoint}{checkpoint object if checkpointing is required.}
  \item{burnin}{burnin size you would like to remove}
  \item{y}{data object}
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
