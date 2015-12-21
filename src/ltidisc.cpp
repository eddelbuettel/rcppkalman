// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// ltidisc.cpp -- based on lti_disc.m

// %LTI_DISC  Discretize LTI ODE with Gaussian Noise
// %
// % Syntax:
// %   [A,Q] = lti_disc(F,L,Qc,dt)
// %
// % In:
// %   F  - NxN Feedback matrix
// %   L  - NxL Noise effect matrix        (optional, default identity)
// %   Qc - LxL Diagonal Spectral Density  (optional, default zeros)
// %   dt - Time Step                      (optional, default 1)
// %
// % Out:
// %   A - Transition matrix
// %   Q - Discrete Process Covariance
// %
// % Description:
// %   Discretize LTI ODE with Gaussian Noise. The original
// %   ODE model is in form
// %
// %     dx/dt = F x + L w,  w ~ N(0,Qc)
// %
// %   Result of discretization is the model
// %
// %     x[k] = A x[k-1] + q, q ~ N(0,Q)
// %
// %   Which can be used for integrating the model
// %   exactly over time steps, which are multiples
// %   of dt.

// % History:
// %   11.01.2003  Covariance propagation by matrix fractions
// %   20.11.2002  The first official version.
// %
// % Copyright (C) 2002, 2003 Simo Särkkä
// %
// % $Id: lti_disc.m 111 2007-09-04 12:09:23Z ssarkka $
// %
// % This software is distributed under the GNU General Public
// % Licence (version 2 or later); please refer to the file
// % Licence.txt, included with the software, for details.

// function [A,Q] = lti_disc(F,L,Q,dt)

//   %
//   % Check number of arguments
//   %
//   if nargin < 1
//     error('Too few arguments');
//   end
//   if nargin < 2
//     L = [];
//   end
//   if nargin < 3
//     Q = [];
//   end
//   if nargin < 4
//     dt = [];
//   end

//   if isempty(L)
//     L = eye(size(F,1));
//   end
//   if isempty(Q)
//     Q = zeros(size(F,1),size(F,1));
//   end
//   if isempty(dt)
//     dt = 1;
//   end

//   %
//   % Closed form integration of transition matrix
//   %
//   A = expm(F*dt);

//   %
//   % Closed form integration of covariance
//   % by matrix fraction decomposition
//   %
//   n   = size(F,1);
//   Phi = [F L*Q*L'; zeros(n,n) -F'];
//   AB  = expm(Phi*dt)*[zeros(n,n);eye(n)];
//   Q   = AB(1:n,:)/AB((n+1):(2*n),:);

#include <RcppArmadillo.h>

arma::mat expm(arma::mat m);    // from expmMat.cpp

//' Discretize Linear Time-Invariant ODE with Gaussian Noise
//'
//' This function discretizes the linear time-invariant (LTI)
//' ordinary differential equation (ODE).
//' @title Discretize Linear Time-Invariant ODE
//' @param F An N x N feedback matrix
//' @param L (Optional, default idendity) N x L noise effect matrix
//' @param Q (Optionalm default zeros) L x L diagonal spectral density
//' @param dt (Option, default one) time step
//' @return A list with elements
//' \describe{
//'   \item{A}{the transition matrix, and}
//'   \item{Q}{the discrete process covariance}
//' }
//' @seealso The documentation for the EKF/UKF toolbox at
//' \url{http://becs.aalto.fi/en/research/bayes/ekfukf}
//' @author The EKF/UKF Toolbox was written by Simo Särkkä, Jouni Hartikainen,
//' and Arno Solin.
//'
//' Dirk Eddelbuettel is porting this package to R and C++, and maintaing it.
// [[Rcpp::export]]
Rcpp::List ltiDisc(const arma::mat & F,
                   const arma::mat & L,
                   const arma::mat & Q,
                   const double dt) {

    // A = expm(F*dt);
    arma::mat A = expm(F * dt);

    // n = size(F,1);
    int n = F.n_rows;

    // Phi = [F L*Q*L'; zeros(n,n) -F'];
    arma::mat Zeros = arma::zeros<arma::mat>(n, n);
    arma::mat Phi = arma::join_cols( arma::join_rows( F,     L * Q * L.t() ),
                                     arma::join_rows( Zeros, -F.t()        ) );

    // AB  = expm(Phi*dt)*[zeros(n,n);eye(n)];
    arma::mat Eye   = arma::eye<arma::mat>(n, n);
    arma::mat AB = expm(Phi*dt) * arma::join_cols(Zeros, Eye);

    // Q   = AB(1:n,:)/AB((n+1):(2*n),:);
    int n1 = n - 1;
    arma::mat lhs = AB.rows(0, n1); 	// arrays starts at zero in C++, corresponds to 1:n
    arma::mat rhs = AB.rows(n, 2*n-1);  // arrays starts at zero in C++, corr. to (n+1):(2*n)
    arma::mat newQ = arma::solve(rhs.t(), lhs.t()).t();	// original code uses '/' aka mrdivive()
                                	// solve() does '\' aka mldivide() -- and using
                                        // the transpose is equivalent problem

    return Rcpp::List::create(Rcpp::Named("A") = A,
                              Rcpp::Named("Q") = newQ);
}
