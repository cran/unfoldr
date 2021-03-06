% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spheroid.R
\name{coefficientMatrixSpheres}
\alias{coefficientMatrixSpheres}
\title{Coefficients for Expectation Maximization algorithm}
\usage{
coefficientMatrixSpheres(bin)
}
\arguments{
\item{bin}{non-decreasing vector of class limits}
}
\value{
array of coefficients
}
\description{
Matrix of coefficients for Wicksell's corpuscle problem
}
\details{
The function computes the matrix of coefficients for solving the discretized integral equation of
the Wicksell's corpuscle problem.
}
\references{
Ohser, J. and Muecklich, F. Statistical analysis of microstructures in materials science J. Wiley & Sons, 2000
}
\author{
M. Baaske
}
