## rcppkalman

Kalman filtering via RcppArmadillo -- based on a port of the 
[EKF/UKF](http://becs.aalto.fi/en/research/bayes/ekfukf/) toolbox for Matlab

### Demos

#### Continuous Wiener-Process Acceleration

See Section 2.2.4 on pages 11 to 15 of the 
[EKF/UKF Documentation](https://github.com/eddelbuettel/rcppkalman/blob/master/inst/ekfukf/ekfukf-documentation.pdf); 
the animation is part of the [corresponding Matlab code](https://github.com/eddelbuettel/rcppkalman/blob/master/inst/ekfukf/demos/kf_cwpa_demo/kf_cwpa_demo.m)

![Smoothing](https://github.com/eddelbuettel/rcppkalman/blob/master/inst/animation/cwpa_smooth.gif)
![Filtering](https://github.com/eddelbuettel/rcppkalman/blob/master/inst/animation/cwpa_filter.gif)

### Status

Still incomplete.  We currently support two basic demo scripts.  Additional functions should get added over time.

### Author

The [EKF/UKF Toolbox for Matlab](http://becs.aalto.fi/en/research/bayes/ekfukf) 
was written by Simo Särkkä, Jouni Hartikainen, and Arno Solin.

Dirk Eddelbuettel is writing and maintaing this packageby porting it to R and
C++ via [Rcpp](https://github.com/RcppCore/Rcpp) and particularly [RcppArmadillo](https://github.com/RcppCore/RcppArmadillo).

### License

This package as well EKF/UKF are released under GNU General Public License,
Version 2 or later. 


