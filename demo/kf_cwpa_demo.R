## % Demonstration for Kalman filter and smoother using a 2D CWPA model
## %
## % Copyright (C) 2007 Jouni Hartikainen
## %
## % This software is distributed under the GNU General Public 
## % Licence (version 2 or later); please refer to the file 
## % Licence.txt, included with the software, for details.

suppressMessages(library(RcppKalman))

kf_cwpa_demo <- function(seed=666) {

    set.seed(seed)
    
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
    Qc <- diag(c(q, q))
    
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
        Xr[, i] <- A %*% Xr[, i-1] + MASS::mvrnorm(1, rep(0, nrow(F)), Q)
    }

    ## % Generate the measurements.
    ## for i = 1:n
    ##    Y(:,i) = H*X_r(:,i) + gauss_rnd(zeros(size(Y,1),1), R);
    ## end
    for (i in 1:n) {
        Y[,i] <- H %*% Xr[,i] + MASS::mvrnorm(1, rep(0, nrow(Y)), R)
    }

    ## clf; clc;
    ## disp('This is a demonstration program for the classical Kalman filter.');
    ## disp(' ');
    ## disp(['KF is used here to estimate the position of a moving object, whos ',...
    ##      'dynamics follow the CWPA-model described in the documentation ',...
    ##       'provided with the toolbox.']);
    ## disp(' ');
    ## disp(['We get noisy measurements from the objects position and velocity ',...
    ##       'on discrete time steps. The real position of the object and the ',...
    ##       'measurements made of them are displayed now. The blue line is the ',...
    ##       'real path of the object and the green dots represents the ',...
    ##       'measurements. The red circle represents the starting point ',...
    ##       'of the object.']);
    ##
    ## disp(' ');
    ## fprintf('Filtering with KF...');
    cat("This is a demonstration program for the classical Kalman filter.\n\n")
    cat("KF is used here to estimate the position of a moving object, whose\n") 
    cat("dynamics follow the CWPA-model described in the documentation\n")
    cat("provided with the toolbox.\n")
    
    ## plot(X_r(1,:),X_r(2,:),Y(1,:),Y(2,:),'.',X_r(1,1),...
    ##      X_r(2,1),'ro','MarkerSize',12);
    ## legend('Real trajectory', 'Measurements');
    ## title('Position');
    plot(Xr[1,], Xr[2,], type='l', main="Position", xlab="", ylab="", col="blue")
    points(Y[1,], Y[2,], col="green", pch=18)
    points(Xr[1,1], Xr[2,1], col='red', pch=1, cex=1.7)
    legend("topleft", c("Real Trajectory", "Measurements"),
           lty=c(1,NA), pch=c(NA,18), col=c("blue", "green"), bty="n")
    
}

kf_cwpa_demo()
