function [var_eps] = AOloopRW(G,H,C_phi_zero,sigmae,phi_sim)

[ sigma, sko ] = AOloop_nocontrol(phi_sim,sigmae,H,G);

usedPhiSim = cell2mat(phi_sim);
[n, m] = size(usedPhiSim);

[U,S,V]     = svd(G,'econ');
r           = rank(G);
S1          = S(1:r,1:r);

% reconstructing wavefront
phik_hat = pinv(G)*sko;

[mG, nG] = size(G);
C = -H;
epsk = zeros(n,m);
epsk_mean_removed = zeros(n,m);
sigma = zeros(m,1);
for i = 1:m-1
    if i == 1
        d = -(C_phi_zero^-1 + G'*G)^-1 * G'*sko(:,1);
        phik_hat(:,i+1) = S1*sko(:,1);
        epsk(:,i+1) = phik_hat(:,i+1);
        epsk_mean_removed(:,i+1) = phik_hat(:,i) - mean(phik_hat(:,i));
    else
        sk = G*(phik_hat(:,i+1) - H*deltaUK);
        epsk(:,i+1) = (C_phi_zero^-1 + G'*G)^-1 * G'*sk - H*deltaUK;
        d = -(C_phi_zero^-1 + G'*G)^-1 * G'*sk;
        phik_hat(:,i+1) = S1*sk;
        epsk_mean_removed(:,i+1) = phik_hat(:,i+1) - mean(phik_hat(:,i+1));
    end
    deltaUK = lsqlin(C,d);
    sigma(i+1) = var(epsk_mean_removed(:,i+1));
end
var_eps = mean(sigma);