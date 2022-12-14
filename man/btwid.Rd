% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decon.R
\name{btwid}
\alias{btwid}
\title{Scaled and Translated B-Spline Fourier Transform}
\usage{
btwid(t, k1, h)
}
\arguments{
\item{t}{Argument of Fourier Transform}

\item{k1}{Leftmost side of basis function}

\item{h}{Knot spacing}
}
\value{
(complex-valued) Fourier transform of scaled and translated B-Spline
}
\description{
Scaled B-Spline Fourier Transform
}
\details{
Evaluates Fourier transform of an integer-knot B-Spline scaled by
h and translated so that the left-most support begins at k1.
}
