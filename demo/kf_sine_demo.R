## % Kalman Filter demonstration with sine signal.
## %
## % History:
## %    3.12.2002 SS  The first implementation
## %
## % Copyright (C) 2002 Simo S"arkk"a
## %
## % This software is distributed under the GNU General Public 
## % Licence (version 2 or later); please refer to the file 
## % Licence.txt, included with the software, for details.

suppressMessages(library(RcppKalman))
#set.seed(42)

  ## %
  ## % Create sine function
  ## %
  ## S1 = [0.2;1.0];
  ## S2 = [1.0;-0.2];
  ## sd = 0.1;
  ## dt = 0.1;
  ## w = 1;
  ## T = (0:dt:30);
  ## X = sin(w*T);
  ## Y = X + sd*randn(size(X));
S1 <- c(0.2, 1.0)
S2 <- c(1.0, -0.2)
stdev <- 0.1
dt <- 0.1
w <- 1
Tseq <- seq(0, 30, by=dt)
X <- sin(w*Tseq)
Y <- X + stdev * rnorm(length(X))

  ## %
  ## % Initialize KF to values
  ## %
  ## %   x = 0
  ## %   dx/dt = 0
  ## %
  ## % with great uncertainty in derivative
  ## %
  ## M = [0;0];
  ## P = diag([0.1 2]);
  ## R = sd^2;
  ## H = [1 0];
  ## q = 0.1;
  ## F = [0 1;
  ##      0 0];
  ## [A,Q] = lti_disc(F,[],diag([0 q]),dt);
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
  ## % Track and animate
  ## %
  ## MM = zeros(size(M,1),size(Y,2));
  ## PP = zeros(size(M,1),size(M,1),size(Y,2));
  ## clf;
  ## clc;
  ## disp('In this demonstration we estimate a stationary sine signal from noisy measurements by using the classical Kalman filter.');
  ## disp(' ');
  ## disp('The filtering results are now displayed sequantially for 10 time step at a time.');
  ## disp(' ');
  ## disp('<push any key to proceed to next time steps>');

n <- length(M)
p <- length(Y)
MM = matrix(0, n, p)
Plist <- vector(length=p, mode="list")
for (i in 1:p) Plist[[i]] <- matrix(0, n, n)


  ## for k=1:size(Y,2)
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
    Plist[[k]] = P

    ## %
    ## % Animate
    ## %
    ## if rem(k,10)==1
    ##   plot(T,X,'b--',...
    ##        T,Y,'ro',...
    ##        T(k),M(1),'k*',...
    ##        T(1:k),MM(1,1:k),'k-');
    ##   legend('Real signal','Measurements','Latest estimate','Filtered estimate')
    ##   title('Estimating a noisy sine signal with Kalman filter.');
    ##   drawnow;
      
    ##   pause;
    ## end
  ## end
}

op <- par(mfcol=c(1,2), mar=c(3,3,1,1))
plot(Tseq, X, type='l', lty="dashed", col="blue", ylim=range(Y))
points(Tseq, Y, col="red", pch="+")
lines(Tseq[1:k], MM[1, 1:k], col="black")
legend("topright", bty="n", lty=c("dashed", NA, "solid"), pch=c(NA, "+", NA),
       legend=c("Real signal", "Measurement", "Filtered estimate"),
       col=c("blue", "red", "black"))


  ## clc;
  ## disp('In this demonstration we estimate a stationary sine signal from noisy measurements by using the classical Kalman filter.');
  ## disp(' ');
  ## disp('The filtering results are now displayed sequantially for 10 time step at a time.');
  ## disp(' ');
  ## disp('<push any key to see the filtered and smoothed results together>')
  ## pause;  
  ## %
  ## % Apply Kalman smoother
  ## %

    
  ## SM = rts_smooth(MM,PP,A,Q);
rl <- rtsSmoother(MM, Plist, A, Q)
SM <- rl[["SM"]]

  ## plot(T,X,'b--',...
  ##      T,MM(1,:),'k-',...
  ##      T,SM(1,:),'r-');
  ## legend('Real signal','Filtered estimate','Smoothed estimate') 
  ## title('Filtered and smoothed estimate of the original signal');

plot(Tseq, X, type='l', lty="dashed", col="blue", ylim=range(Y))
lines(Tseq, MM[1, ], col="black")
lines(Tseq, SM[1, ], col="red")
legend("topright", bty="n", lty=c("dashed", "solid", "solid"), 
       legend=c("Real signal", "Filtered", "Smoothed"),
       col=c("blue", "black", "red"))
par(op)


  ## clc;
  ## disp('The filtered and smoothed estimates of the signal are now displayed.')
  ## disp(' ');
  ## disp('RMS errors:');
  ## %
  ## % Errors
  ## %
  ## fprintf('KF = %.3f\nRTS = %.3f\n',...
  ##         sqrt(mean((MM(1,:)-X(1,:)).^2)),...
  ##         sqrt(mean((SM(1,:)-X(1,:)).^2)));

cat("RMS errors for KF and RTS:\n")
print(sqrt(mean( (MM[1,] - X)^2)))
print(sqrt(mean( (SM[1,] - X)^2)))

