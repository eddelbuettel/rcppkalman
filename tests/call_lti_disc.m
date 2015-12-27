# -*- mode: Matlab; -*- 

## usage: call_lti_disc()
function A = call_lti_disc()

    addpath("../inst/ekfukf")
    addpath("../RcppKalman/ekfukf")     

    process_noise_variance = 0.01;

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
    A
    Q

endfunction
