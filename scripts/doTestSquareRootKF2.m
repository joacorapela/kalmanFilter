
loadRes = load('~/dev/research/programs/github/python/pykalman/pykalman/datasets/data/robot.mat');
% Initialize simulation variables
SRSigmaW = chol(loadRes.Q,'lower'); % Square-root process noise covar
SRSigmaV = chol(loadRes.R,'lower'); % Square-root sensor noise covar
A = loadRes.A; B = zeros(size(A,1), 1); C = loadRes.C; D = zeros(size(C,1), 1); % Plant definition matrices
N = size(loadRes.x, 2);
xHat0 = loadRes.x0; % Initialize estimated system initial state
SigmaX0 = loadRes.P_0; % Initialize Kalman filter initial covariance
% u = 0; % Unknown initial driving input: assume zero
us = zeros([1, N]);
ws = SRSigmaW*randn([size(A, 1), N]);
vs = SRSigmaV*randn([size(C, 1), N]);
xs = zeros(size(A, 1), N);
zs = zeros(size(C, 1), N);
for k=2:N
    xs(:,k)=A*xs(:,k-1) + B*us(:,k) + ws(:,k);
    zs(:,k)=C*xs(:,k) + D*us(:,k) + vs(:,k);
end
[xHats, SigmaXHats] = squareRootKF(A, B, C, D, xHat0, SigmaX0, SRSigmaW, SRSigmaV, us, zs, N);

% figure(1); clf;
% plot(0:N-1,xs(1:N),'k-',0:N-1,xHats,'b--', ...
% 0:N-1,xHats+3*sqrt(SigmaXHats),'m-.',...
% 0:N-1,xHats-3*sqrt(SigmaXHats),'m-.'); grid;
% title('Kalman filter in action'); xlabel('Iteration');
% ylabel('State'); legend('true','estimate','bounds');
% figure(2); clf;
% plot(0:N-1,xs(1:N)-xHats,'b-', ...
% 0:N-1,3*sqrt(SigmaXHats),'m--',...
% 0:N-1,-3*sqrt(SigmaXHats),'m--');
% % grid; legend('Error','bounds',0); title('Error with bounds');
% grid; legend('Error','bounds'); title('Error with bounds');
% xlabel('Iteration'); ylabel('Estimation Error');
