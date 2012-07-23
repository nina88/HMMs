\name{equil}
\alias{hmm.sim}
\alias{equil}
\alias{rdiric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
equil
}
\description{
calculates equilibrium distribution, random number generator for dirichlet distribution, simulates hidden markov model
}
\usage{
hmm.sim(n, lambda, P, f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
length of observed sequence
}
  \item{lambda}{
hidden sequence transition matrix
}
  \item{P}{
array of transition matrices for observed sequence
}
  \item{f}{
number of states of observed sequence
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{s }{segmentation}
%%  \item{y }{observed sequence}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nina Wilkinson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n, lambda, P, f) 
{
    r = dim(lambda)[1]
    pi.lam = equil(lambda)
    s = numeric(n)
    s[1] = sample(1:r, 1, replace = TRUE, prob = pi.lam)
    for (t in 2:n) {
        s[t] = sample(1:r, 1, replace = TRUE, prob = lambda[match(s[t - 
            1], 1:r), ])
    }
    pi.P = array(0, c(f, r))
    for (i in 1:r) {
        pi.P[, i] = equil(P[, , i])
    }
    y = numeric(n)
    y[1] <- sample(1:f, 1, replace = TRUE, prob = pi.P[, s[1]])
    for (t in 2:n) {
        y[t] <- sample(1:f, 1, replace = TRUE, prob = P[y[t - 
            1], , s[t]])
    }
    list(s = s, y = y)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
