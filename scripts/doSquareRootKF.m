% Initialize simulation variables
SRSigmaW = chol(1,'lower'); % Square-root process noise covar
SRSigmaV = chol(1,'lower'); % Square-root sensor noise covar
A = 1; B = 1; C = 1; D = 0; % Plant definition matrices
maxIter = 40;
xtrue = 0; xhat = 0 % Initialize true and estimated system initial state
SigmaX = 0.1; % Initialize Kalman filter covariance
SRSigmaX = chol(SigmaX,'lower');
u = 0; % Unknown initial driving input: assume zero
% Reserve storage for variables we might want to plot/evaluate
xstore = zeros(maxIter+1,length(xtrue)); xstore(1,:) = xtrue;
xhatstore = zeros(maxIter,length(xhat));
SigmaXstore = zeros(maxIter,length(xhat)); % store diagonal only
for k = 1:maxIter,
% SR-KF Step 1a: State estimate time update
xhat = A*xhat + B*u; % use prior value of "u"
% SR-KF Step 1b: Error covariance time update
SRSigmaX = qr([A*SRSigmaX, SRSigmaW]')';
SRSigmaX = tril(SRSigmaX(1:length(xhat),1:length(xhat)));
% [Implied operation of system in background, with
% input signal u, and output signal z]
u = 0.5*randn(1) + cos(k/pi); % for example... (measured)
w = SRSigmaW*randn(length(xtrue));
v = SRSigmaV*randn(length(C*xtrue));
ztrue = C*xtrue + D*u + v; % y is based on present x and u
xtrue = A*xtrue + B*u + w; % future x is based on present u
% SR-KF Step 1c: Estimate system output
zhat = C*xhat + D*u;
% SR-KF Step 2a: Compute Kalman gain matrix
% Note: "help mrdivide" to see how "division" is implemented
SRSigmaZ = qr([C*SRSigmaX,SRSigmaV]')';
SRSigmaZ = tril(SRSigmaZ(1:length(zhat),1:length(zhat)));
L = (SRSigmaX*SRSigmaX')*C'/SRSigmaZ'/SRSigmaZ;
% SR-KF Step 2b: State estimate measurement update
xhat = xhat + L*(ztrue - zhat);
% SR-KF Step 2c: Error covariance measurement update
Sx_ = SRSigmaX';
cov_update_vectors = L*SRSigmaZ;
for j=1:length(zhat),
Sx_ = cholupdate(Sx_,cov_update_vectors(:,j),'-');
end
SRSigmaX = Sx_';
% [Store information for evaluation/plotting purposes]
xstore(k+1,:) = xtrue; xhatstore(k,:) = xhat;
SigmaXstore(k,:) = diag(SRSigmaX*SRSigmaX');
end;
figure(1); clf;
plot(0:maxIter-1,xstore(1:maxIter),'k-',0:maxIter-1,xhatstore,'b--', ...
0:maxIter-1,xhatstore+3*sqrt(SigmaXstore),'m-.',...
0:maxIter-1,xhatstore-3*sqrt(SigmaXstore),'m-.'); grid;
title('Kalman filter in action'); xlabel('Iteration');
ylabel('State'); legend('true','estimate','bounds');
figure(2); clf;
plot(0:maxIter-1,xstore(1:maxIter)-xhatstore,'b-', ...
0:maxIter-1,3*sqrt(SigmaXstore),'m--',...
0:maxIter-1,-3*sqrt(SigmaXstore),'m--');
% grid; legend('Error','bounds',0); title('Error with bounds');
grid; legend('Error','bounds'); title('Error with bounds');
xlabel('Iteration'); ylabel('Estimation Error');
