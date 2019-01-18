function [ A, Cw, K ] = computeKalmanAR(C_phi_zero, C_phi_one, G, sigmae)
    % computing A
    A = C_phi_one / C_phi_zero;
    computedC_phi1 = A*C_phi_zero;
    % computing Cw
    Cw = (eye-A^2)*C_phi_zero;
    R = eye(72);
    Q = (Cw+Cw.')/2;
    disp(size(G'))
    P = dare(A',G',Q,R);
    % compute Kalman gain
    K = (S + A*P*G')*(G*P*G' + sigmae*eye)^-1;
end