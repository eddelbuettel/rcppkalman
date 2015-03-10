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
