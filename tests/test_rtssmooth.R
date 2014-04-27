
suppressMessages(library(xts))
suppressMessages(library(RcppKalman))
suppressMessages(library(RcppOctave))

setwd("~/git/rcppkalman/tests/")

SPY <- readRDS("spy.rds")
colnames(SPY) <- c("Open", "High", "Low", "Close", "Volume", "Adjusted")

o_source(file = "call_rts_smooth.m")
mP <- .CallOctave("call_rts_smooth", SPY)
#cat("mP after Octave is\n")

## copied / adapted from funPart1
call_rts_smooth <- function(data) {

    #tick <- 0.25            # with tick <- 0.25 results identical ... but should be 0.01
    #tick <- 0.025

    #n = length(data(:,4));
    n <- nrow(data)
    
    ##finish = input('enter finish, no greater than n  ')
    ##if ( finish > length(data(:,4)) )
    ##   finish = 0 % i.e. all available data is used
    ##end
    finish <- 0

    ##open = data(:,4) ;
    ##high = data(:,5) ;
    ##low = data(:,6) ;
    ##close = data(:,7) ;
    ##market_type = data(:,230) ;
    open <- data[,1]
    high <- data[,2]
    low  <- data[,3]
    close <- data[,4]
    market_type <- 1  #; #data(:,230) ;

    #clear data

    #vwap = round( ( ( open .+ close .+ ( (high .+ low) ./ 2 ) ) ./ 3 ) ./ tick) .* tick ;
    #vwap <- round( ( ( open + close + ( (high + low) / 2 ) ) / 3 ) / tick) * tick
    vwap <- ( open + close + ( (high + low) / 2 ) ) / 3 
    
    #vwap_process_noise = ( vwap .- shift(vwap,1) ) ./ 2.0 ;
    #vwap_process_noise = ( tail(vwap, -1) - head(vwap, -1) ) / 2.0 

    #median_vwap_process_noise = median(vwap_process_noise(2:end-finish,1)) ;
    #median_vwap_process_noise <- median(vwap_process_noise[-1])
    
    #vwap_process_noise_deviations = vwap_process_noise(2:end-finish,1) .- median_vwap_process_noise ;
    #vwap_process_noise_deviations <- vwap_process_noise[-1] - median_vwap_process_noise 
    
    #MAD_process_noise <- median( abs( vwap_process_noise_deviations ) ) 

    ## convert this to variance under the assumption of a normal distribution
    #std_vwap_noise <- 1.4826 * MAD_process_noise ;
    #process_noise_variance <- std_vwap_noise * std_vwap_noise ;

    #measurement_noise <- 0.666 * ( high - low ) ;
    #median_measurement_noise <- median( measurement_noise ) #(1:end-finish,1) ) ;
    ##measurement_noise_deviations = measurement_noise(1:end-finish,1) .- median_measurement_noise ;
    #measurement_noise_deviations <- measurement_noise - median_measurement_noise 
    #MAD_measurement_noise <- median( abs( measurement_noise_deviations ) ) 

    ## convert this to variance under the assumption of a normal distribution
    #std_measurement_noise <- 1.4826 * MAD_measurement_noise ;
    #measurement_noise_variance <- std_measurement_noise * std_measurement_noise ;

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
    #Qc = diag([q q]);
    Qc <- diag(2) * process_noise_variance 
    #print(Qc)
    
    ## Discretisation of the continuous-time system.
    #[A,Q] = lti_disc(F,L,Qc,1); % last item is dt stepsize set to 1

    #print(F)
    #print(L)
    #print(Qc)
    rl <- ltiDisc(F, L, Qc, 1)
    A <- rl[["A"]]
    Q <- rl[["Q"]]
    #print(Q)
    #A


    ## Measurement model.
    ##H = [1 0 0 0 0 0;
    ##     0 1 0 0 0 0];
    H <- matrix(0, 2, 6)
    H[1,1] = H[2,2] = 1
    
    ## Variance in the measurements.
    #r1 = measurement_noise_variance ;
    #R = diag([r1 r1]);
    r1 <- measurement_noise_variance
    R <- diag(2) * measurement_noise_variance
    #print(R)
    
    ## Initial guesses for the state mean and covariance
    m <- matrix(c(0, vwap[1], 0, 0, 0, 0), 6, 1)
    P <- diag(c(0.1, 0.1, 0.1, 0.1, 0.5, 0.5))
              
    ## Space for the estimates.
    #MM = zeros(size(m,1), length(vwap));
    n <- length(vwap)              
    MM <- matrix(0, 6, n)
    #MM <- vector(length=n, mode="list")

    ## create vectors for eventual plotting
    #predict_plot = zeros(length(vwap),1) ;
    #MM_plot = zeros(length(vwap),1) ;
    #sigmaP_plus = zeros(length(vwap),1) ;
    #sigmaP_minus = zeros(length(vwap),1) ;

    ## see help for rts_smooth
    #PP = zeros(size(m,1), size(m,1), length(vwap));
    PP <- vector(length=n, mode="list")              

    B <- diag(6)
    u <- matrix(0, 6, 1)
    
    #for ii = 1:length(vwap)
    for (ii in 1:n) {

        #[m,P] = kf_predict(m,P,A,Q);
        #if (ii == 1) print(P)
        #m <- kfPredict(m, P, A, Q, B, u)
        rl <- kfPredict(m, P, A, Q, B, u)
        m <- rl[["x"]]
        P <- rl[["P"]]
        #if (ii == 1825) { cat("after kfPred\n"); print(P) }
        ## #print(m)

        #[m,P] = kf_update(m,P,vwap(ii,:),H,R);
        #m <- kfUpdate(m, P, vwap[ii], H, R)
        rl <- kfUpdate(m, P, vwap[ii], H, R)
        m <- rl[["x"]]
        P <- rl[["P"]]
        #if (ii == 1825) { cat("after kfUpd\n"); print(P) }

        #if (ii < 4) {
        #    print(P);
        #}
        K <- rl[["K"]]
        
        ## if (ii <= 3) {
        ## #    print(m)
        ## #    stop("now")
        ##     cat("After kfUpdate, P is now\n")
        ##     print(P)
        ## }
        
        ## store MM and PP here
        MM[, ii] <- m
        #PP[[ii]] <- P
        ## there is potential R bug here and we force a copy in this function
        ##PP <- storeInList(PP, P, ii-1);
        PP[[ii]] <- P
        #if (ii == 1825) { # || ii == 5) {
        ##     print(class(PP[[ii]]))
        #    cat("PP_i")
        #    print(PP[[ii]])
        #}
        
        ## sigmaP is for storing the current error covariance for plotting purposes
        #sigmaP <- sqrt(diag(P)) ; 
        #sigmaP_plus(ii,1) = MM_plot(ii,1) + 2 * sigmaP(1) ;
        #sigmaP_minus(ii,1) = MM_plot(ii,1) - 2 * sigmaP(1) ;
        
    }
    #print(m)
    #cat("In R, MM and P index 1\n")
    #print(MM[,1,drop=FALSE])
    #print(PP[[1]])
    #[SM,SP] = rts_smooth(MM,PP,A,Q);
    ## cat("After loop first three P are\n")
    ## print(class(PP))
    ## print(class(PP[[1]]))
    ## #print(str(PP))
    ## for (i in 1:3) {
    ##     print(PP[[i]])
    ## }

    #for (i in 1834:1836) print(PP[[i]])

    #print(MM[,1,drop=FALSE])
    #cat("After loop\n")
    #print(PP[[1825]])
    #print(MM[,5,drop=FALSE])
    #print(PP[[5]])
    
    rl <- rtsSmoother(MM, PP, A, Q)
}


rl <- call_rts_smooth(coredata(SPY))
SM <- rl[["SM"]]
rP <- SM[, ncol(SM), drop=FALSE]
SP <- rl[["SP"]]
print(rP)
print(SM[,1831:1836])
print(diag(SP[[1834]]))
print(diag(SP[[1835]]))
print(diag(SP[[1836]]))
print( all.equal(mP,rP))



