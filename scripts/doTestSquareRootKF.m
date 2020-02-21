
% Initialize simulation variables
SRSigmaW = chol(1,'lower'); % Square-root process noise covar
SRSigmaV = chol(1,'lower'); % Square-root sensor noise covar
A = 1; B = 1; C = 1; D = 0; % Plant definition matrices
maxIter = 40;
xhat0 = 0; % Initialize estimated system initial state
SigmaX0 = 0.1; % Initialize Kalman filter initial covariance
% u = 0; % Unknown initial driving input: assume zero
us = 0.5*randn([1, maxIter])+cos((1:maxIter)/pi);
ws = SRSigmaW*randn(length(us));
vs = SRSigmaV*randn(length(C*us));
xs = zeros(maxIter, 1);
zs = zeros(maxIter, 1);
for k=2:maxIter
    xs(k)=A*xs(k-1) + B*us(k) + ws(k);
    zs(k)=C*xs(k) + D*us(k) + vs(k);
end
[xhats, SigmaXHats] = squareRootKF(A, B, C, D, xhat0, SigmaX0, SRSigmaW, SRSigmaV, us, zs, maxIter);

figure(1); clf;
plot(0:maxIter-1,xs(1:maxIter),'k-',0:maxIter-1,xhats,'b--', ...
0:maxIter-1,xhats+3*sqrt(SigmaXHats),'m-.',...
0:maxIter-1,xhats-3*sqrt(SigmaXHats),'m-.'); grid;
title('Kalman filter in action'); xlabel('Iteration');
ylabel('State'); legend('true','estimate','bounds');
figure(2); clf;
plot(0:maxIter-1,xs(1:maxIter)-xhats,'b-', ...
0:maxIter-1,3*sqrt(SigmaXHats),'m--',...
0:maxIter-1,-3*sqrt(SigmaXHats),'m--');
% grid; legend('Error','bounds',0); title('Error with bounds');
grid; legend('Error','bounds'); title('Error with bounds');
xlabel('Iteration'); ylabel('Estimation Error');
