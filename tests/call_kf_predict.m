# -*- mode: Matlab; -*- 

## usage: call_kf_predict()
function m = call_kf_predict()

    addpath("../inst/ekfukf")

    process_noise_variance = 0.01;
    measurement_noise_variance = 0.01;

    ## Transition matrix for the continous-time system.
    F = [0 0 1 0 0 0;
         0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1;
         0 0 0 0 0 0;
         0 0 0 0 0 0];

    ## Noise effect matrix for the continous-time system.
    L = [0 0;
         0 0;
         0 0;
         0 0;
         1 0;
         0 1];

    ## Process noise variance
    q = process_noise_variance ;
    Qc = diag([q q]);

    ## Discretisation of the continuous-time system.
    [A,Q] = lti_disc(F,L,Qc,1); % last item is dt stepsize set to 1

    ## Measurement model.
    H = [1 0 0 0 0 0;
         0 1 0 0 0 0];

    ## Variance in the measurements.
    r1 = measurement_noise_variance ;
    R = diag([r1 r1]);

    ## Initial guesses for the state mean and covariance
    #m = [0 vwap(1,1) 0 0 0 0]';
    m = [0 150.75 0 0 0 0]';
    P = diag([0.1 0.1 0.1 0.1 0.5 0.5]) ;

    [m,P] = kf_predict(m,P,A,Q);
    m
    P
    
endfunction
