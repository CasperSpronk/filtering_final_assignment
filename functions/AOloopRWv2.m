function [var_eps] = AOloopRW(G,H,C_phi_zero,sigmae,phi_sim)

[ sigma, sko ] = AOloop_nocontrol(phi_sim,sigmae,H,G);

usedPhiSim = cell2mat(phi_sim);
dataPoints = length(usedPhiSim);

[U,S,V]     = svd(G,'econ');
r           = rank(G);
U1          = U(:,1:r);
S1          = S(1:r,1:r);
V1          = V(:,1:r);

% reconstructing wavefront
phik = V1*S*U1*sko;

phikplus1 = zeros(dataPoints);
for i = 2:dataPoints
    phikplus1(i) = phik(i-1) + wgn(1,n);
end




