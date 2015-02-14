## % Demonstration for Kalman filter and smoother using a 2D CWPA model
## %
## % Copyright (C) 2007 Jouni Hartikainen
## %
## % This software is distributed under the GNU General Public 
## % Licence (version 2 or later); please refer to the file 
## % Licence.txt, included with the software, for details.

suppressMessages(library(RcppKalman))

kf_cwpa_demo <- function() {

    ## % Transition matrix for the continous-time system.
    ## F = [0 0 1 0 0 0;
    ##      0 0 0 1 0 0;
    ##      0 0 0 0 1 0;
    ##      0 0 0 0 0 1;
    ##      0 0 0 0 0 0;
    ##      0 0 0 0 0 0];
    F <- matrix(0, 6, 6)
    F[1,3] <- F[2,4] <- F[3,5] <- F[4,6] <- 1
        
    ## % Noise effect matrix for the continous-time system.
    ## L = [0 0;
    ##      0 0;
    ##      0 0;
    ##      0 0;
    ##      1 0;
    ##      0 1];
    L <- matrix(0, 6, 2)
    L[5,1] <- L[6,2] <- 1
        
    ## % Stepsize
    ## dt = 0.5;
    dt <- 0.5
    
    ## % Process noise variance
    ## q = 0.2;
    ## Qc = diag([q q]);
    q <- 0.2
    QC <- diag(c(q, q))
    
    ## % Discretization of the continous-time system.
    ## [A,Q] = lti_disc(F,L,Qc,dt);
    rl <- ltiDisc(F, L, Qc, dt)
    A <- rl[["A"]]
    Q <- rl[["Q"]]

    ## % Measurement model.
    ## H = [1 0 0 0 0 0;
    ##      0 1 0 0 0 0];
    H <- matrix(0, 2, 6)
    H[1,1] <- H[2,2] <- 1
    
    ## % Variance in the measurements.
    ## r1 = 10;
    ## r2 = 5;
    ## R = diag([r1 r1]);
    r1 <- 10
    r2 <- 5
    R <- diag(c(r1, r2))

    ## % Generate the data.
    ## n = 50;
    ## Y = zeros(size(H,1),n);
    ## X_r = zeros(size(F,1),n);
    ## X_r(:,1) = [0 0 0 0 0 0]';
    ## for i = 2:n
    ##    X_r(:,i) = A*X_r(:,i-1) + gauss_rnd(zeros(size(F,1),1), Q);
    ## end
    n <- 50
    Y <- matrix(0, nrow(H), n)
    Xr <- matrix(0, nrow(F), n)
    for (i in 2:n) {
        Xr[, i] <- A %*% Xr[, i-1] + MASS::mvrnorn(1, rep(0, nrow(F)), Q)
    }
    
    
}
