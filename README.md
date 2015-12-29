## rcppkalman

Kalman filtering via RcppArmadillo -- based on a R and C++ port of the 
[EKF/UKF](http://becs.aalto.fi/en/research/bayes/ekfukf/) toolbox for Matlab

### Demos

#### Static Sine Signal And Noisy Measurement

This example is not described in the pdf manual, but included as demo
[kf_sine_demo.m](https://github.com/eddelbuettel/rcppkalman/blob/master/inst/ekfukf/demos/kf_sine_demo/kf_sine_demo.m) within the
[EKF/UKF](http://becs.aalto.fi/en/research/bayes/ekfukf/) sources. A signal
is provided via a sine wave plus random noise, and a linear Kalman Filter is
used to smooth and filter the series. Our variant
[demo/kf_sine_demo.R](https://github.com/eddelbuettel/rcppkalman/blob/master/demo/kf_sine_demo.R)
reproduces the demo via the following chart

![Sine Signal](https://github.com/eddelbuettel/rcppkalman/blob/master/inst/images/kf_sine_demo.png)

#### Continuous Wiener-Process Acceleration

This demo is described in detail in Section 2.2.4 on pages 11 to 15 of the 
[EKF/UKF Documentation](https://github.com/eddelbuettel/rcppkalman/blob/master/inst/ekfukf/ekfukf-documentation.pdf); 
the animation is part of the corresponding Matlab code in
[kf_cwpa_demo.m](https://github.com/eddelbuettel/rcppkalman/blob/master/inst/ekfukf/demos/kf_cwpa_demo/kf_cwpa_demo.m). We
show the two final charts which provide animations of the smoothing and
filtering in our version [demo/kf_cwpa_demo.R](https://github.com/eddelbuettel/rcppkalman/blob/master/demo/kf_cwpa_demo.R):

![Smoothing](https://github.com/eddelbuettel/rcppkalman/blob/master/inst/animation/cwpa_smooth.gif)
![Filtering](https://github.com/eddelbuettel/rcppkalman/blob/master/inst/animation/cwpa_filter.gif)

### Status

Working, but still far from complete.  We currently support two demo scripts
based on linear smoothers and filters.  Additional functions should get added
over time.

### Author

The [EKF/UKF Toolbox for Matlab](http://becs.aalto.fi/en/research/bayes/ekfukf) 
was written by Simo Särkkä, Jouni Hartikainen, and Arno Solin.

Dirk Eddelbuettel is writing and maintaing this package by porting it to R and
C++ via [Rcpp](https://github.com/RcppCore/Rcpp) and particularly [RcppArmadillo](https://github.com/RcppCore/RcppArmadillo).

### License

This package is released under GNU General Public License, Version 2 or
later. EKF/UKF itself (which is included) is released under the GNU General
Public License, Versions 2 and 3.


