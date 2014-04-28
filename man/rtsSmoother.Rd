\name{rtsSmoother}
\alias{rtsSmoother}
\title{Rauch-Tung-Striebel smoother}
\usage{
rtsSmoother(M, Plist, A, Q)
}
\arguments{
  \item{M}{An N x K matrix of K mean estimates from the
  Kalman Filter}

  \item{Plist}{A list of length N with N x N matrices of
  state covariances from the Kalman Filter}

  \item{A}{An N x N state transition matrix (or in the more
  general case a list of K such matrices; not yet
  implemented)}

  \item{Q}{An N x N noise covariance matrix (or in the more
  general case a list of K such matrices; not yet
  implemented)}
}
\value{
A list with three elements M, P and D which are,
respectively, the smoothed mean sequence, the smooted state
covariance sequence, and the smoothed gain sequence.
}
\description{
The function computes the Rauch-Tung-Striebel smoother.
}
\details{
This function implements the Rauch-Tung-Striebel smoother
algorithm which calculate a "smoothed" sequence from the
given Kalman filter output sequence by conditioning all
steps to all measurements.
}
\author{
The EKF/UKF Toolbox was written by Simo Särkkä, Jouni
Hartikainen, and Arno Solin. Dirk Eddelbuettel is writing
and maintaing this package by porting it to R and C++.
}
\seealso{
The EKF/UKF toolbox at
\url{http://becs.aalto.fi/en/research/bayes/ekfukf}
}

