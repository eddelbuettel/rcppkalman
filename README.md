## rcppkalman

Kalman filtering via RcppArmadillo -- based on a R and C++ port of the 
[EKF/UKF](http://becs.aalto.fi/en/research/bayes/ekfukf/) toolbox for Matlab

[![Build Status](https://travis-ci.org/eddelbuettel/rcppkalman.svg)](https://travis-ci.org/eddelbuettel/rcppkalman) 
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html) 
[![Last Commit](https://img.shields.io/github/last-commit/eddelbuettel/rcppkalman)](https://github.com/eddelbuettel/rcppkalman)

## Why?

R has a number of existing Kalman filter packages which are all very
noteworthy in their own right. A comparison was provided by
[this JSS paper from 2011](http://www.jstatsoft.org/v39/i02).

Yet I had a need for something both simple and fast at the C++ level.

The [EKF/UKF](http://becs.aalto.fi/en/research/bayes/ekfukf/) toolbox for
Matlab proved to be a wonderful source of excellent code that was well
documented (see [this 130 page pdf manual]http://becs.aalto.fi/en/research/bayes/ekfukf/documentation.pdf()),
under a suitable license and covering both simple examples as
well as promising extensions.

### Demos

#### Static Sine Signal And Noisy Measurement

This example is not described in the pdf manual, but included as demo
[kf_sine_demo.m](https://github.com/eddelbuettel/rcppkalman/blob/master/inst/ekfukf/demos/kf_sine_demo/kf_sine_demo.m) within the
[EKF/UKF](http://becs.aalto.fi/en/research/bayes/ekfukf/) sources. A signal
is provided via a sine wave plus random noise, and a linear Kalman Filter is
used to smooth and filter the series. Our variant
[demo/kf_sine_demo.R](https://github.com/eddelbuettel/rcppkalman/blob/master/demo/kf_sine_demo.R)
reproduces the demo via the following chart

![Sine Signal](http://eddelbuettel.github.io/rcppkalman/kf_sine_demo.png)

#### Continuous Wiener-Process Acceleration

This demo is described in detail in Section 2.2.4 on pages 11 to 15 of the 
[EKF/UKF Documentation](https://github.com/eddelbuettel/rcppkalman/blob/master/inst/ekfukf/ekfukf-documentation.pdf); 
the animation is part of the corresponding Matlab code in
[kf_cwpa_demo.m](https://github.com/eddelbuettel/rcppkalman/blob/master/inst/ekfukf/demos/kf_cwpa_demo/kf_cwpa_demo.m). We
show the two final charts which provide animations of the smoothing and
filtering in our version [demo/kf_cwpa_demo.R](https://github.com/eddelbuettel/rcppkalman/blob/master/demo/kf_cwpa_demo.R):

![Smoothing](http://eddelbuettel.github.io/rcppkalman/cwpa_smooth.gif)
![Filtering](http://eddelbuettel.github.io/rcppkalman/cwpa_filter.gif)

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


