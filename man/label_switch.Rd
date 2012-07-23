\name{label_switch}
\alias{label_switch}
\alias{permut}
\alias{order_of_firsts}
\alias{permutations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
label switching
}
\description{
label switching
} 
\usage{
label_switch(segment.store2, segment.store1, r, lambda, P, f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segment.store2}{
segmentation at iteration i
}
  \item{segment.store1}{
segmentation at iteration i-1
}
  \item{r}{
number of segments
}
  \item{lambda}{
transition matrices for segments
}
  \item{P}{
array of transition matrices for observed sequence
}
  \item{f}{
number of states in observed sequence
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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
function (segment.store2, segment.store1, r, lambda, P, f) 
{
    x = numeric(factorial(r))
    new.P = array(0, c(f, f, r))
    new.lambda = matrix(0, nrow = r, ncol = r)
    for (i in 1:factorial(r)) {
        perms = permutation(r, segment.store2, i)
        x[i] = sum(perms == segment.store1)
    }
    y = which.max(x)
    segment.store3 = permutation(r, segment.store2, y)
    z1 = order.of.firsts(segment.store2, r)
    z2 = order.of.firsts(segment.store3, r)
    if (identical(segment.store2, segment.store3)) {
        new.P = P
        new.lambda = lambda
    }
    else {
        for (k in 1:r) {
            new.P[, , z2[k]] = P[, , z1[k]]
            for (l in 1:r) {
                new.lambda[z2[k], z2[l]] = lambda[z1[k], z1[l]]
            }
        }
    }
    return(list(segment.store2 = segment.store3, P = new.P, lambda = new.lambda))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
