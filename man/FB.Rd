\name{FB}
\alias{FB}
\title{The forward--backward algorithm}
\description{
The forward--backward algorithm
}
\usage{
FB(y, lambda, P)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{hmm_fasta object}
  \item{lambda}{transition matrix for hidden sequence}
  \item{P}{array of transition matrices for each hidden state}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{s }{segmentation}
  \item{xi }{normalising constant}
  \item{back}{backwards probabilities}
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
function (y, lambda, P) 
{
    n = length(y)
    r = dim(lambda)[1]
    h = dim(p)[1]
    f = matrix(0, nrow = r, ncol = n)
    s = numeric(n)
    back = array(0, c(r, r, n))
    xi = numeric(n)
    lml = 0
    pi.P = array(0, c(h, r))
    for (i in 1:r) {
        pi.P[, i] = equil(p[, , i])
    }
    pi.lambda = equil(lambda)
    xi[1] = sum(pi.P[y[1], ] * pi.lambda)
    f[, 1] = (pi.P[y[1], ] * pi.lambda)/xi[1]
    lml = lml + log(xi[1])
    for (i in 2:n) {
        for (k in 1:r) {
            f[k, i] = p[y[i - 1], y[i], k] * sum(lambda[, k] * 
                f[, i - 1])
        }
        xi[i] = sum(f[, i])
        f[, i] = f[, i]/xi[i]
        lml = lml + log(xi[i])
    }
    for (i in (n - 1):1) {
        back[, , i] = lambda * f[, i]
        for (k in 1:r) {
            back[, k, i] = back[, k, i]/sum(back[, k, i])
        }
    }
    s[n] = sample(1:r, 1, prob = f[, n])
    for (i in (n - 1):1) {
        s[i] = sample(1:r, 1, prob = back[, s[i + 1], i])
    }
    return(list(s = s, xi = xi, back = back))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
