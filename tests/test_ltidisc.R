
suppressMessages(library(xts))
suppressMessages(library(RcppKalman))
if (requireNamespace("RcppOctave", quietly=TRUE)) {

    suppressMessages(library(RcppOctave))

    o_source(file = "call_lti_disc.m")
    mA <- .CallOctave("call_lti_disc")

    ## copied / adapted from funPart1.m
    call_lti_disc <- function() {

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
        print(A)
        print(Q)
        A
    }


    rA <- call_lti_disc()

    print( all.equal(mA,rA))
}
