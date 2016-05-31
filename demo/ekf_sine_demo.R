
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

ekf_sine_demo <- function() {

}
