// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// rtssmooth.cpp -- based on rts_smooth.m

// %RTS_SMOOTH  Rauch-Tung-Striebel smoother
// %
// % Syntax:
// %   [M,P,S] = RTS_SMOOTH(M,P,A,Q)
// %
// % In:
// %   M - NxK matrix of K mean estimates from Kalman filter
// %   P - NxNxK matrix of K state covariances from Kalman Filter
// %   A - NxN state transition matrix or NxNxK matrix of K state
// %       transition matrices for each step.
// %   Q - NxN process noise covariance matrix or NxNxK matrix
// %       of K state process noise covariance matrices for each step.
// %
// % Out:
// %   M - Smoothed state mean sequence
// %   P - Smoothed state covariance sequence
// %   D - Smoother gain sequence
// %
// % Description:
// %   Rauch-Tung-Striebel smoother algorithm. Calculate "smoothed"
// %   sequence from given Kalman filter output sequence
// %   by conditioning all steps to all measurements.
// %
// % Example:
// %   m = m0;
// %   P = P0;
// %   MM = zeros(size(m,1),size(Y,2));
// %   PP = zeros(size(m,1),size(m,1),size(Y,2));
// %   for k=1:size(Y,2)
// %     [m,P] = kf_predict(m,P,A,Q);
// %     [m,P] = kf_update(m,P,Y(:,k),H,R);
// %     MM(:,k) = m;
// %     PP(:,:,k) = P;
// %   end
// %   [SM,SP] = rts_smooth(MM,PP,A,Q);
// %
// % See also:
// %   KF_PREDICT, KF_UPDATE

// % Copyright (C) 2003-2006 Simo S?rkk?
// %
// % $Id: rts_smooth.m 109 2007-09-04 08:32:58Z jmjharti $
// %
// % This software is distributed under the GNU General Public
// % Licence (version 2 or later); please refer to the file
// % Licence.txt, included with the software, for details.

// function [M,P,D] = rts_smooth(M,P,A,Q)

//   %
//   % Check which arguments are there
//   %
//   if nargin < 4
//     error('Too few arguments');
//   end

//   %
//   % Extend A and Q if they are NxN matrices
//   %
//   if size(A,3)==1
//     A = repmat(A,[1 1 size(M,2)]);
//   end
//   if size(Q,3)==1
//     Q = repmat(Q,[1 1 size(M,2)]);
//   end

//   %
//   % Run the smoother
//   %
//   D = zeros(size(M,1),size(M,1),size(M,2));
//   for k=(size(M,2)-1):-1:1
//     P_pred   = A(:,:,k) * P(:,:,k) * A(:,:,k)' + Q(:,:,k);
//     D(:,:,k) = P(:,:,k) * A(:,:,k)' / P_pred;
//     M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - A(:,:,k) * M(:,k));
//     P(:,:,k) = P(:,:,k) + D(:,:,k) * (P(:,:,k+1) - P_pred) * D(:,:,k)';
//   end

#include <RcppArmadillo.h>



//' The function computes the Rauch-Tung-Striebel smoother.
//'
//' This function implements the Rauch-Tung-Striebel smoother algorithm which
//' calculate a "smoothed" sequence from the given Kalman filter output sequence
//' by conditioning all steps to all measurements.
//'
//' @title Rauch-Tung-Striebel smoother
//' @param M An N x K matrix of K mean estimates from the Kalman Filter
//' @param P An N x N x K cube length K with N x N state covariances matrices 
//' from the Kalman Filter
//' @param A An N x N state transition matrix (or in the more general case a
//' list of K such matrices; not yet implemented)
//' @param Q An N x N noise covariance matrix  (or in the more general case a
//' list of K such matrices; not yet implemented)
//' @return A list with three elements
//' \describe{
//'   \item{SM}{the smoothed mean sequence,}
//'   \item{SP}{the smooted state covariance sequence,and }
//'   \item{D}{the smoothed gain sequence.}
//' }
//' @author The EKF/UKF Toolbox was written by Simo Särkkä, Jouni Hartikainen,
//' and Arno Solin.
//'
//' Dirk Eddelbuettel is porting this package to R and C++, and maintaing it.
//' @seealso The documentation for the EKF/UKF toolbox at
//' \url{http://becs.aalto.fi/en/research/bayes/ekfukf}
// [[Rcpp::export]]
Rcpp::List rtsSmoother(arma::mat & M, 		
                       arma::cube & P,          
                       const arma::mat & A,
                       const arma::mat & Q) {
    
    int n = M.n_rows;
    int k = M.n_cols;

    //   %
    //   % Extend A and Q if they are NxN matrices
    //   %
    //   if size(A,3)==1
    //     A = repmat(A,[1 1 size(M,2)]);
    //   end
    //   if size(Q,3)==1
    //     Q = repmat(Q,[1 1 size(M,2)]);
    //   end
    arma::cube AA(n,n,k), QQ(n,n,k);
    for (int i=0; i<k; i++) {
        AA.slice(i) = A;
        QQ.slice(i) = Q;
    }


    //   %
    //   % Run the smoother
    //   %
    //   D = zeros(size(M,1),size(M,1),size(M,2));
    //   for k=(size(M,2)-1):-1:1
    //     P_pred   = A(:,:,k) * P(:,:,k) * A(:,:,k)' + Q(:,:,k);
    //     D(:,:,k) = P(:,:,k) * A(:,:,k)' / P_pred;
    //     M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - A(:,:,k) * M(:,k));
    //     P(:,:,k) = P(:,:,k) + D(:,:,k) * (P(:,:,k+1) - P_pred) * D(:,:,k)';
    //   end
    arma::cube D = arma::zeros<arma::cube>(n, n, k);

    for (unsigned int j=k-1-1; j>0; j--) {
        arma::mat Ppred = AA.slice(j) * P.slice(j) * AA.slice(j).t() + QQ.slice(j);
        arma::mat lhs = (P.slice(j) * AA.slice(j).t());
        arma::mat rhs = Ppred;
        D.slice(j) = arma::solve(rhs.t(), lhs.t()).t();
        M.col(j) = M.col(j) + D.slice(j) * (M.col(j+1) - AA.slice(j) * M.col(j));
        P.slice(j) = P.slice(j) + D.slice(j) * (P.slice(j+1) - Ppred) * D.slice(j).t();
    }

    return Rcpp::List::create(Rcpp::Named("SM") = M,
                              Rcpp::Named("SP") = P,
                              Rcpp::Named("D") = D);

}




