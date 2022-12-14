% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createBasis.R
\name{createBasis}
\alias{createBasis}
\title{createBasis}
\usage{
createBasis(rangeval, qdim)
}
\arguments{
\item{rangeval}{Numeric vector of length 2 indicating the interval
[a, b] on which the basis may be evaluated.}

\item{qdim}{Integer dimension of basis}
}
\value{
A list, where `$B`, `$Bp`, and `$Bpp` evaluate the basis
    functions and their first two derivatives, respectively, `$G`
    and `$P` are inner-product matrices for `$B` and `$Bpp`,
    respectively, and `knots` and `h` are the knots and
    knot-spacing associated with the basis.
}
\description{
Constructs q-dimensional basis on an interval [a, b]
}
\author{
David Kent
}
