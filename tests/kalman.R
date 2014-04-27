
## cf test_kfpredict.R
##
kalman <- function(data, tick=0.01) {

    n <- nrow(data)
    
    finish <- 0

    open <- data[,1]
    high <- data[,2]
    low  <- data[,3]
    close <- data[,4]
    market_type <- 1 

    vwap <- round( ( ( open + close + ( (high + low) / 2 ) ) / 3 ) / tick) * tick
    ## vwap_process_noise = ( tail(vwap, -1) - head(vwap, -1) ) / 2.0 
    ## median_vwap_process_noise <- median(vwap_process_noise[-1])
    ## vwap_process_noise_deviations <- vwap_process_noise[-1] - median_vwap_process_noise 
    ## MAD_process_noise <- median( abs( vwap_process_noise_deviations ) ) 

    ## ## convert this to variance under the assumption of a normal distribution
    ## std_vwap_noise <- 1.4826 * MAD_process_noise 
    ## process_noise_variance <- std_vwap_noise * std_vwap_noise 

    ## measurement_noise <- 0.666 * ( high - low ) 
    ## median_measurement_noise <- median( measurement_noise ) #(1:end-finish,1) ) ;
    ## measurement_noise_deviations <- measurement_noise - median_measurement_noise 
    ## MAD_measurement_noise <- median( abs( measurement_noise_deviations ) ) 

    ## ## convert this to variance under the assumption of a normal distribution
    ## std_measurement_noise <- 1.4826 * MAD_measurement_noise 
    ## measurement_noise_variance <- std_measurement_noise * std_measurement_noise 

    process_noise_variance <- 0.01
    measurement_noise_variance <- 0.01

    
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
    q <- process_noise_variance 
    Qc <- diag(2) * process_noise_variance 

    ## Discretisation of the continuous-time system.
    #[A,Q] = lti_disc(F,L,Qc,1); % last item is dt stepsize set to 1
    Q <- matrix(0, 6, 6);
    A <- ltiDisc(F, L, Qc, 1, Q)

    ## Measurement model.
    ##H = [1 0 0 0 0 0;
    ##     0 1 0 0 0 0];
    H <- matrix(0, 2, 6)
    H[1,1] = H[2,2] = 1
    
    ## Variance in the measurements.
    r1 <- measurement_noise_variance
    R <- diag(2) * measurement_noise_variance
    
    ## Initial guesses for the state mean and covariance
    m <- matrix(c(0, vwap[1], 0, 0, 0, 0), 6, 1)
    P <- diag(c(0.1, 0.1, 0.1, 0.1, 0.5, 0.5))
              
    ## Space for the estimates.
    MM <- matrix(0, 6, n)

    ## create vectors for eventual plotting
    #predict_plot = zeros(length(vwap),1) ;
    #MM_plot = zeros(length(vwap),1) ;
    #sigmaP_plus = zeros(length(vwap),1) ;
    #sigmaP_minus = zeros(length(vwap),1) ;

    ## see help for rts_smooth
    PP <- vector(length=n, mode="list")              

    B <- diag(6)
    u <- matrix(0, 6, 1)
    
    for (ii in 1:n) {

        #[m,P] = kf_predict(m,P,A,Q);
        m <- kfPredict(m, P, A, Q, B, u)

        #[m,P] = kf_update(m,P,vwap(ii,:),H,R);
        m <- kfUpdate(m, P, vwap[ii], H, R)
        
        ## store MM and PP here
        MM[, ii] <- m
        PP[[ii]] <- P
        
        ## sigmaP is for storing the current error covariance for plotting purposes
        #sigmaP <- sqrt(diag(P)) ; 
        #sigmaP_plus(ii,1) = MM_plot(ii,1) + 2 * sigmaP(1) ;
        #sigmaP_minus(ii,1) = MM_plot(ii,1) - 2 * sigmaP(1) ;
    }
    
    #[SM,SP] = rts_smooth(MM,PP,A,Q);
    SM <- rtsSmoother(MM, PP, A, Q)
    
    invisible(list(m=m, SM=SM))
}
