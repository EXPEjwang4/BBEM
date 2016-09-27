\name{BBEM_bs}
\alias{BBEM_bs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian bootstrap EM algorithm for variable selection in linear model
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BBEM_bs(y, X, n.var, v0 = 0.1, v1 = 100, maxIter = 20, T = 100, threshold = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  A vector of response variable.
%%     ~~Describe \code{y} here~~
}
  \item{X}{
  Design matrix which has the same number of rows as the length of y.
%%     ~~Describe \code{X} here~~
}
  \item{n.var}{
  Number of variables that are used in each boostrap sample. Should not exceed the number of columns in the design matrix.
%%     ~~Describe \code{n.var} here~~
}
  \item{v0}{
  Hyperparameter in the spike and slab prior for the coefficient. Default is 0.1.
%%     ~~Describe \code{v0} here~~
}
  \item{v1}{
  Hyperparameter in the spike and slab prior for the coefficient. Default is 100.
%%     ~~Describe \code{v1} here~~
}
  \item{maxIter}{
  Maximum number of iterations. Default is 20.
%%     ~~Describe \code{maxIter} here~~
}
  \item{T}{
  Number of Bootstrap samples. Default is 100.
%%     ~~Describe \code{T} here~~
}
  \item{threshold}{
  Number of iterations that the estimated model index stays the same before the algorithm stops. Default is 3.
%%     ~~Describe \code{threshold} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{gamma}{The average estimation of model index \eqn{\gamma} from T bootstrap samples}
  \item{theta}{The estimated parameter \eqn{\theta} in the Binomial distribution}
  \item{sigmaSq}{The estimated variance parameter \eqn{\sigma^2}}
  \item{selvar}{The number of times each variable being selected}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
http://arxiv.org/abs/1603.04360
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
