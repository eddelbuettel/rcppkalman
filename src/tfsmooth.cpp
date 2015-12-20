// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// tfsmooth.cpp -- based on tf_smooth.m

// %TF_SMOOTH  Two filter based Smoother
// %
// % Syntax:
// %   [M,P] = TF_SMOOTH(M,P,Y,A,Q,H,R,[use_inf])
// %
// % In:
// %   M - NxK matrix of K mean estimates from Kalman filter
// %   P - NxNxK matrix of K state covariances from Kalman Filter
// %   Y - Sequence of K measurement as DxK matrix
// %   A - NxN state transition matrix.
// %   Q - NxN process noise covariance matrix.
// %   H - DxN Measurement matrix.
// %   R - DxD Measurement noise covariance.
// %   use_inf - If information filter should be used (default 1)
// %
// % Out:
// %   M - Smoothed state mean sequence
// %   P - Smoothed state covariance sequence
// %   
// % Description:
// %   Two filter linear smoother algorithm. Calculate "smoothed"
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
// %   [SM,SP] = tf_smooth(MM,PP,A,Q,H,R,Y);
// %
// % See also:
// %   KF_PREDICT, KF_UPDATE

// % History:
// %   
// %   02.8.2007 JH Changed the name to tf_smooth
// %   26.3.2007 JH Fixed a bug in backward filter with observations having
// %                having more than one dimension.
// %             
// % Copyright (C) 2006 Simo Särkkä
// %               2007 Jouni Hartikainen
// %
// % $Id: tf_smooth.m 111 2007-09-04 12:09:23Z ssarkka $
// %
// % This software is distributed under the GNU General Public 
// % Licence (version 2 or later); please refer to the file 
// % Licence.txt, included with the software, for details.
// %

// function [M,P] = tf_smooth(M,P,Y,A,Q,H,R,use_inf)

#include <RcppArmadillo.h>

Rcpp::List kfUpdate(const arma::vec & x, const arma::mat & P, const arma::vec & y,
                    const arma::mat & H, const arma::mat & R);
Rcpp::List kfPredict(const arma::vec & x, const arma::mat & P, const arma::mat & A,
                     const arma::mat & Q, const arma::mat & B, const arma::vec & u);

// [[Rcpp::export]]
Rcpp::List tfSmoother(const arma::mat & Mc, 		
                      const arma::cube & Pc,          
                      const arma::mat & Y,
                      const arma::mat & A,
                      const arma::mat & Q,
                      const arma::mat & H,
                      const arma::mat & R,
                      const bool useinf) {

    arma::mat M = Mc;
    arma::cube P = Pc;
    
    int n = M.n_rows;
    int k = M.n_cols;
    // %
    // % Run the backward filter
    // %
    // if use_inf
    //   zz = zeros(size(M));
    //   SS = zeros(size(P));
    //   IR = inv(R);
    //   IQ = inv(Q);
    //   z = zeros(size(M,1),1);
    //   S = zeros(size(M,1),size(M,1));
    //   for k=size(M,2):-1:1
    //     G = S / (S + IQ);
    //     S = A' * (eye(size(M,1)) - G) * S * A;
    //     z = A' * (eye(size(M,1)) - G) * z;
    //     zz(:,k)   = z;
    //     SS(:,:,k) = S;
    //     S = S + H'*IR*H;
    //     z = z + H'*IR*Y(:,k);
    //   end
    // else
    //   BM = zeros(size(M));
    //   BP = zeros(size(P));
    //   IA = inv(A);
    //   IQ = IA*Q*IA';  
    //   fm = zeros(size(M,1),1);
    //   fP = 1e12*eye(size(M,1));
    //   BM(:,end) = fm;
    //   BP(:,:,end) = fP;
    //   for k=(size(M,2)-1):-1:1
    //     [fm,fP] = kf_update(fm,fP,Y(:,k+1),H,R);
    //     [fm,fP] = kf_predict(fm,fP,IA,IQ);
    //     BM(:,k) = fm;
    //     BP(:,:,k) = fP;
    //   end
    // end
    arma::mat zz = arma::zeros<arma::mat>(n, k);
    arma::cube SS = arma::zeros<arma::cube>(n, n, k);
    arma::mat BM = arma::zeros<arma::mat>(n, k);
    arma::cube BP = arma::zeros<arma::cube>(n, n, k);
    if (useinf) {
        arma::mat IR = arma::inv(R);
        arma::mat IQ = arma::inv(Q);
        arma::vec z = arma::zeros<arma::vec>(n);
        arma::mat S = arma::zeros<arma::mat>(n, n);
        for (int j=k-1; j>=0; j--) {
            arma::mat rhs = S + IQ;
            arma::mat lhs = S;
            arma::mat G = arma::solve(rhs.t(), lhs.t()).t();  // cf ltidisc.cpp for discussion of solve
            S = A.t() * (arma::eye(n,n) - G) * S * A;
            z = A.t() * (arma::eye(n,n) - G) * z;
            zz.col(j) = z;
            SS.slice(j) = S;
            S = S + H.t() * IR * H;
            z = z + H.t() * IR * Y.col(j);
        }
    } else {
        arma::mat IA = arma::inv(A);
        arma::mat IQ = IA * Q * IA.t();
        arma::vec fm = arma::zeros<arma::vec>(n);
        arma::vec fP = 1e12 * arma::eye(n,n);
        BM.col(k) = fm;
        BP.slice(k) = fP;
        arma::mat B = arma::eye(n,n); // to add to kfPredict -- FIXME
        arma::vec u = arma::zeros<arma::vec>(n);
        for (int j=k-2; j>=0; j--) {
            Rcpp::List rl = kfUpdate(fm, fP, Y.col(j), H, R);
            fm = Rcpp::as<arma::vec>(rl[0]);
            fP = Rcpp::as<arma::vec>(rl[1]);
            rl = kfPredict(fm, fP, IA, IQ, B, u);
            fm = Rcpp::as<arma::vec>(rl[0]);
            fP = Rcpp::as<arma::vec>(rl[1]);
            BM.col(j) = fm;
            BP.slice(j) = fP;
        }
    }
    
    // %
    // % Combine estimates
    // %
    // if use_inf
    //   for k=1:size(M,2)-1
    //     G = P(:,:,k) * SS(:,:,k) / (eye(size(M,1)) + P(:,:,k) * SS(:,:,k));
    //     P(:,:,k) = inv(inv(P(:,:,k)) + SS(:,:,k));
    //     M(:,k) = M(:,k) + P(:,:,k) * zz(:,k) - G * M(:,k);
    //   end
    // else
    //   for k=1:size(M,2)-1
    //     tmp = inv(inv(P(:,:,k)) + inv(BP(:,:,k)));
    //     M(:,k) = tmp * (P(:,:,k)\M(:,k) + BP(:,:,k)\BM(:,k));
    //     P(:,:,k) = tmp;
    //   end
    // end
    if (useinf) {
        for (int j=0; j<k-1; j++) {
            arma::mat rhs = arma::eye(n, n) + P.slice(j) * SS.slice(j);
            arma::mat lhs = P.slice(j) * SS.slice(j);
            arma::mat G = arma::solve(rhs.t(), lhs.t()).t();  // cf ltidisc.cpp for discussion of solve
            P.slice(j) = arma::inv(arma::inv(P.slice(j)) + SS.slice(j));
            M.col(j) = M.col(j) + P.slice(j) * zz.col(j) - G * M.col(j);
        }
    } else {
        for (int j=0; j<k-1; j++) {
            arma::mat tmp = arma::inv(arma::inv(P.slice(j) + inv(BP.slice(j))));
            M.col(j) = tmp * (arma::solve(P.slice(j), M.col(j)) + arma::solve(BP.slice(j), BM.col(j)));
            P.slice(j) = tmp;
        }
    }

    return Rcpp::List::create(Rcpp::Named("M") = M,
                              Rcpp::Named("P") = P);
    
}
