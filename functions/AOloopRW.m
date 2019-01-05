function [var_eps] = AOloopRW(G,H,C_phi_zero,sigma_e,phi_sim)

[ sigma, sk ] = AOloop_nocontrol(phi_sim,sigmae,H,G)

usedPhiSim = cell2mat(phi_sim);
% Taking Singular Value Decomposition of G
[U,S,V]     = svd(G,'econ');                                                   
% Partitioning SVD as in (Verhaegen, p.34)
r           = rank(G);
U1          = U(:,1:r);
S1          = S(1:r,1:r);
V1          = V(:,1:r);
% Expressing estimate of the wavefront
phi_est     = V1*inv(S1)*U1'*sk;