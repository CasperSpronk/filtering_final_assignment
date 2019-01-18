function [var_eps] = AOloopRWwithS(G,H,C_phi_zero,sigmae,phi_sim)

[ sigma, sko ] = AOloop_nocontrol(phi_sim,sigmae,H,G);

usedPhiSim = cell2mat(phi_sim);
[n, m] = size(usedPhiSim);

% Taking partitioned svd of G
[U,S,V]     = svd(G,'econ');
r           = rank(G);
U1          = U(:,1:r);
S1          = S(1:r,1:r);
V1          = V(:,1:r);

% reconstructing wavefront
phik_hat = V1*S1^(-1)*U1'*sko;

u = H^-1*phik_hat;
epskplus1 = zeros(size(usedPhiSim));
epsk_meanless = zeros(size(usedPhiSim));
for i = 1:length(usedPhiSim)-1
    epskplus1(:,i+1) = phik_hat(:,i+1) - H*u(:,i);
    epsk_meanless(:,i+1) = epskplus1(:,i+1) - mean(epskplus1(:,i+1));
    var_epsk_meanless = var(epsk_meanless(:,i+1));
end

var_eps = mean(var_epsk_meanless);





