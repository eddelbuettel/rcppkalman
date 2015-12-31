
suppressMessages(library(RcppKalman))
if ((Sys.info()[["sysname"]] != "Windows") &&	 ## win-builder has non-working RcppOctave, so exclude
    (requireNamespace("RcppOctave", quietly=TRUE))) {

    suppressMessages(library(RcppOctave))

    #setwd("~/git/rcppkalman/tests/")

    ## data generation as in the demo 'kf_sine_demo'
    set.seed(42)
    S1 <- c(0.2, 1.0)
    S2 <- c(1.0, -0.2)
    stdev <- 0.1
    dt <- 0.1
    w <- 1
    Tfinal <- 30
    Tseq <- seq(0, Tfinal, by=dt)
    X <- sin(w*Tseq)
    Y <- X + stdev * rnorm(length(X))

    o_source(file = "call_tf_smooth.m")
    MM <- .CallOctave("call_tf_smooth", Y)
    mP <- MM[, ncol(MM), drop=FALSE]

    call_tf_smooth <- function(Y) {

        ##print(Y[1:10])
        ## %
        ## % Initialize KF to values
        ## %
        ## %   x = 0
        ## %   dx/dt = 0
        ## %
        ## % with great uncertainty in derivative
        ## %
        M <- c(0,0)
        P <- diag(c(0.1, 2))
        R <- matrix(stdev^2,1,1)
        H <- matrix(c(1, 0), 1, 2)
        q <- 0.1
        Qc <- diag(c(0, q))
        F <- matrix(0, 2, 2)
        F[1,2] <- 1
        L <- diag(2)
        rl <- ltiDisc(F, L, Qc, dt)
        A <- rl[["A"]]
        Q <- rl[["Q"]]

        ## %
        ## % Track 
        ## %
        ## MM = zeros(size(M,1),size(Y,2));
        ## PP = zeros(size(M,1),size(M,1),size(Y,2));
        n <- length(M)
        p <- length(Y)
        MM <- matrix(0, n, p)
        PP <- array(0, dim=c(n, n, p))

        for (k in 1:p) {
            
            ## %
            ## % Track with KF
            ## %
            ## [M,P] = kf_predict(M,P,A,Q);
            ## [M,P] = kf_update(M,P,Y(k),H,R);
            B <- diag(n)
            u <- matrix(0, n, 1)
            rl <- kfPredict(M, P, A, Q, B, u)
            M <- rl[["x"]]
            P <- rl[["P"]]

            rl <- kfUpdate(M, P, Y[k], H, R)
            M <- rl[["x"]]
            P <- rl[["P"]]

            MM[,k] <- M
            PP[,,k] <- P
        }
        rl <- rtsSmoother(MM, PP, A, Q)
        rl <- tfSmoother(MM, PP, matrix(Y,nrow=1), A, Q, H, R, TRUE)
        rl
    }

    rl <- call_tf_smooth(Y)
    SM <- rl[["M"]]
    rP <- SM[, ncol(SM), drop=FALSE]
    print( all.equal(mP,rP))
    print( all.equal(MM, SM)) 

}
