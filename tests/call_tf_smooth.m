# -*- mode: Matlab; -*- 

## usage: call_tf_smooth(data)
function SM2 = call_tf_smooth(Y)

    ## also needs some Kalman Filter code -- from http://becs.aalto.fi/en/research/bayes/ekfukf/
    addpath("../inst/ekfukf")
    addpath("../RcppKalman/ekfukf")     

    n = length(Y(:));

    dt = 0.1;                           %  both from data generation
    sd = 0.1;
    
    %
    % Initialize KF to values
    %
    %   x = 0
    %   dx/dt = 0
    %
    % with great uncertainty in derivative
    %
    M = [0;0];
    P = diag([0.1 2]);
    R = sd^2;
    H = [1 0];
    q = 0.1;
    F = [0 1;
         0 0];
    [A,Q] = lti_disc(F,[],diag([0 q]),dt);

    MM = zeros(size(M,1),size(Y,2));
    PP = zeros(size(M,1),size(M,1),size(Y,2));

    for k=1:size(Y,2)
        %
        % Track with KF
        %
        [M,P] = kf_predict(M,P,A,Q);
        [M,P] = kf_update(M,P,Y(k),H,R);
    
        MM(:,k) = M;
        PP(:,:,k) = P;
    end
    %
    % Apply Kalman smoother
    %
    % Smoothing step.
    [SM,SP] = rts_smooth(MM,PP,A,Q);
    [SM2,SP2] = tf_smooth(MM,PP,Y,A,Q,H,R,1);
    %MM(:,297:301)
    %SM(:,297:301)
    
endfunction


