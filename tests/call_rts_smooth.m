# -*- mode: Matlab; -*- 

## usage: call_rts_smooth(data)
function m = call_rts_smooth(data)

    ## also needs some Kalman Filter code -- from http://becs.aalto.fi/en/research/bayes/ekfukf/
    addpath("../inst/ekfukf")

    n = length(data(:,4));
    finish = 0;

    open = data(:,1) ;
    high = data(:,2) ;
    low = data(:,3) ;
    close = data(:,4) ;

    vwap = ( open .+ close .+ ( (high .+ low) ./ 2 ) ) ./ 3 ;
    # vwap_process_noise = ( vwap .- shift(vwap,1) ) ./ 2.0 ;

    # median_vwap_process_noise = median(vwap_process_noise(2:end-finish,1)) ;
    # vwap_process_noise_deviations = vwap_process_noise(2:end-finish,1) .- median_vwap_process_noise ;
    # MAD_process_noise = median( abs( vwap_process_noise_deviations ) ) ;

    # ## convert this to variance under the assumption of a normal distribution
    # std_vwap_noise = 1.4826 * MAD_process_noise ;
    # process_noise_variance = std_vwap_noise * std_vwap_noise ;

    # measurement_noise = 0.666 .* ( high .- low ) ;
    # median_measurement_noise = median( measurement_noise(1:end-finish,1) ) ;
    # measurement_noise_deviations = measurement_noise(1:end-finish,1) .- median_measurement_noise ;
    # MAD_measurement_noise = median( abs( measurement_noise_deviations ) ) ;

    # ## convert this to variance under the assumption of a normal distribution
    # std_measurement_noise = 1.4826 * MAD_measurement_noise ;
    # measurement_noise_variance = std_measurement_noise * std_measurement_noise ;

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
    m = [0 vwap(1,1) 0 0 0 0]';
    P = diag([0.1 0.1 0.1 0.1 0.5 0.5]);


    ## Space for the estimates.
    MM = zeros(size(m,1), length(vwap));

    ## create vectors for eventual plotting
    predict_plot = zeros(length(vwap),1) ;
    MM_plot = zeros(length(vwap),1) ;
    sigmaP_plus = zeros(length(vwap),1) ;
    sigmaP_minus = zeros(length(vwap),1) ;

    ## see help for rts_smooth
    PP = zeros(size(m,1), size(m,1), length(vwap));

    for ii = 1:length(vwap)

        [m,P] = kf_predict(m,P,A,Q);

        predict_plot(ii,1) = m(2,1) ;

        [m,P] = kf_update(m,P,vwap(ii,:),H,R);

        MM(:,ii) = m;
        PP(:,:,ii) = P;

        MM_plot(ii,1) = m(2,1) ;

        ## sigmaP is for storing the current error covariance for plotting purposes
        sigmaP = sqrt(diag(P)) ; 
        sigmaP_plus(ii,1) = MM_plot(ii,1) + 2 * sigmaP(1) ;
        sigmaP_minus(ii,1) = MM_plot(ii,1) - 2 * sigmaP(1) ;
    end

    [SM,SP] = rts_smooth(MM,PP,A,Q);
    disp "m, and last four SM and (S)P in Octave"
    SM(:,1833:1836)
    diag(SP(:,:,1833))'
    diag(SP(:,:,1834))'
    diag(SP(:,:,1835))'
    diag(SP(:,:,1836))'

endfunction


