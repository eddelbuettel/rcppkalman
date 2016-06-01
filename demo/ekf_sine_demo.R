
## based on inst/ekfukf/demos/ekf_sine_demo/*

#% Demonstration for EKF using a random sine signal model. 
#%
#% A Very simple demonstration for extended Kalman filter (EKF), which is
#% used to track a random single-component sinusoid signal,
#% which is modelled as x_k = a_k*sin(\theta_k), dtheta/dt = omega_k.
#% The signal is also filtered with unscented Kalman filter (UKF) for
#% comparison.
#%
#% Copyright (C) 2007 Jouni Hartikainen
#%
#% This software is distributed under the GNU General Public 
#% Licence (version 2 or later); please refer to the file 
#% Licence.txt, included with the software, for details.

suppressMessages(library(RcppKalman))

ekf_sine_f <- function(x, param) {
    ##% Dynamical model function for the random sine signal demo

    ##% Copyright (C) 2007 Jouni Hartikainen
    ##%
    ##% This software is distributed under the GNU General Public 
    ##% Licence (version 2 or later); please refer to the file 
    ##% Licence.txt, included with the software, for details.

    ##function x_n = ekf_sine_f(x,param)
    ##dt = param(1);
    ##A = [1 dt 0;0 1 0;0 0 1];
    ##x_n = A*x(1:3,:);
    ##if size(x,1) == 7 || size(x,1) == 6
    ##    x_n(1:3,:) = x_n(1:3,:) + x(4:6,:);
    ##end
    dt <- param[1]
    A <- matrix(c(1, dt, 0,
                  0,  1, 0,
                  0,  0, 1), nrow=3, byrow=TRUE)
    x_n <- A %*% x[1:3, ]
    if (nrow(x) == 7 || nrow(x) == 6) {
        x_n[1:3, ] <- x_n[1:3,] + x <- n[4:6,]
    }
    x_n
}

ekf_sine_h <- function(x, param) {
    ##% Measurement model function for the random sine signal demo

    ##% Copyright (C) 2007 Jouni Hartikainen
    ##%
    ##% This software is distributed under the GNU General Public 
    ##% Licence (version 2 or later); please refer to the file 
    ##% Licence.txt, included with the software, for details.

    ##function Y = ekf_sine_h(x,param)
    ##f = x(1,:);
    ##a = x(3,:);

    ##Y = a.*sin(f);
    ##if size(x,1) == 7
    ##  Y = Y + x(7,:);
    ##end
    f <- x[1, ]
    a <- x[3, ]
    Y <- a * sin(f)
    if (nrow(x) == 7) {
        Y <- Y + x[7,]
    }
    Y
}

ekf_sine_demo <- function() {

    ##clc;
    ##disp('Filtering the signal with EKF...');

    ##save_plots = 1;

    ##% Measurement model and it's derivative
    ##f_func = @ekf_sine_f;
    ##h_func = @ekf_sine_h;
    ##dh_dx_func = @ekf_sine_dh_dx;
    ##d2h_dx2_func = @ekf_sine_d2h_dx2;
    f_func <- ekf_sine_f
    h_func <- ekf_sine_h

    invisible(NULL)
}

