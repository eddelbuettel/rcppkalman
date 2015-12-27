
suppressMessages(library(RcppKalman))
if (requireNamespace("RcppOctave", quietly=TRUE)) {

    suppressMessages(library(RcppOctave))

    o_source(file = "call_kf_update.m")
    mA <- .CallOctave("call_kf_update")

    ## copied / adapted from funPart1.m
    call_kf_update <- function() {

        ## Transition matrix for the continous-time system.
        ## F = [0 0 1 0 0 0;
        ##      0 0 0 1 0 0;
        ##      0 0 0 0 1 0;
        ##      0 0 0 0 0 1;
        ##      0 0 0 0 0 0;
        ##      0 0 0 0 0 0];
        F <- matrix(0, 6, 6)
        F[1,3] <- F[2,4] <- F[3,5] <- F[4,6] <- 1
        
        
        ## Noise effect matrix for the continous-time system.
        ## L = [0 0;
        ##      0 0;
        ##      0 0;
        ##      0 0;
        ##      1 0;
        ##      0 1];
        L <- matrix(0, 6, 2)
        L[5,1] <- L[6,2] <- 1

        ## Process noise variance
        process_noise_variance  <- 0.01 #
                                        #Qc = diag([q q]);
        Qc <- diag(2) * process_noise_variance 

        ## Discretisation of the continuous-time system.
                                        #[A,Q] = lti_disc(F,L,Qc,1); % last item is dt stepsize set to 1

                                        #print(F)
                                        #print(L)
                                        #print(Qc)
        rl <- ltiDisc(F, L, Qc, 1)
        A <- rl[["A"]]
        Q <- rl[["Q"]]

        
        ## Measurement model.
        ##H = [1 0 0 0 0 0;
        ##     0 1 0 0 0 0];
        H <- matrix(0, 2, 6)
        H[1,1] = H[2,2] = 1
        
        ## Variance in the measurements.
                                        #r1 = measurement_noise_variance ;
                                        #R = diag([r1 r1]);
        measurement_noise_variance  <- 0.01 #
        r1 <- measurement_noise_variance
        R <- diag(2) * measurement_noise_variance
        
        ## Initial guesses for the state mean and covariance
                                        #m <- matrix(c(0, vwap[1], 0, 0, 0, 0), 6, 1)
        m <- matrix(c(0, 150.75, 0, 0, 0, 0), 6, 1)
        P <- diag(c(0.1, 0.1, 0.1, 0.1, 0.5, 0.5))
        
        B <- diag(6)
        u <- matrix(0, 6, 1)
        
        ##[m,P] = kf_predict(m,P,A,Q);
        ##if (ii == 1) print(P)
        rl <- kfPredict(m, P, A, Q, B, u)
        m <- rl[["x"]]
        P <- rl[["P"]]

        ##[m,P] = kf_update(m,P,vwap(ii,:),H,R);
        rl <- kfUpdate(m, P, rep(150.75,2), H, R)
        m <- rl[["x"]]
        P <- rl[["P"]]
        K <- rl[["K"]]
        print(m)
        print(P)
                                        #print(K)
        m
    }


    rA <- call_kf_update()

    print( all.equal(mA,rA))
}
