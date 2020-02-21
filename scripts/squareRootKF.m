function [xHats, SigmaXHats] = squareRootKF(A, B, C, D, xHat0, SigmaX0, SRSigmaW, SRSigmaV, us, zs, maxIter)
    xHat = xHat0;
    SRSigmaX = chol(SigmaX0, 'lower');
    xHats = zeros(length(xHat0), maxIter);
    SigmaXHats = zeros(length(xHat), maxIter); % store diagonal only
    for k = 1:maxIter,
        z = zs(:,k);
        u = us(:,k);
        % SR-KF Step 1a: State estimate time update
        xHat = A*xHat + B*u; % use prior value of "u"
        % SR-KF Step 1b: Error covariance time update
        SRSigmaX = qr([A*SRSigmaX, SRSigmaW]')';
        SRSigmaX = tril(SRSigmaX(1:length(xHat),1:length(xHat)));
        % SR-KF Step 1c: Estimate system output
        zhat = C*xHat + D*u;
        % SR-KF Step 2a: Compute Kalman gain matrix
        % Note: "help mrdivide" to see how "division" is implemented
        SRSigmaZ = qr([C*SRSigmaX,SRSigmaV]')';
        SRSigmaZ = tril(SRSigmaZ(1:length(zhat),1:length(zhat)));
        L = (SRSigmaX*SRSigmaX')*C'/SRSigmaZ'/SRSigmaZ;
        % stat debug
        lhs = SRSigmaZ';
        rhs = (SRSigmaX*SRSigmaX')*C';
        % L = rhs/lhs;
        L2 = (lhs'\rhs')';
        L1 = (SRSigmaX*SRSigmaX')*C'/SRSigmaZ';
        % end debug
        % SR-KF Step 2b: State estimate measurement update
        xHat = xHat + L*(z - zhat);
        % SR-KF Step 2c: Error covariance measurement update
        Sx_ = SRSigmaX';
        cov_update_vectors = L*SRSigmaZ;
        for j=1:length(zhat),
            Sx_ = cholupdate(Sx_,cov_update_vectors(:,j),'-');
        end
        SRSigmaX = Sx_';
        % [Store information for evaluation/plotting purposes]
        xHats(:, k) = xHat;
        SigmaXHats(:, k) = diag(SRSigmaX*SRSigmaX');
    end;
