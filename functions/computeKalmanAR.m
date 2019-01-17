function [ A, Cw, K ] = computeKalmanAR(C_phi_zero, C_phi_one, G, sigmae)
    % computing A
    A = C_phi_one / C_phi_zero;
    computedC_phi1 = A*C_phi_zero;
    if computedC_phi1 == C_phi_one
        Cw = (eye-A^2)*C_phi_zero;
        R = sigmae*eye;
    	Q = Cw;
        P = dare(A',G',Q,R);
        K = (S + A*P*G')*(G*P*G' + sigmae*eye)^-1;
        disp("correct")
    else
        disp("wrong")
        K = 0;
        Cw = 0;
    end
    % computing Cw


end