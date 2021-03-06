function [var_eps] = AOloopAR(G,H,C_phi0,sigmae,A,Cw,K,phiSim(1));

[sigmaeNoControl, sko] = AOloop_nocontrol(phiSim(1),sigmae,H,G);
usedPhiSim = cell2mat(phiSim(1));
ukmin1 = 0;
eps_hat = 0;
C = -H;
shatk = sko(:,1);
eps_hat_mean_removed = zeros(size(usedPhiSim));
for i = 2:length(usedPhiSim)
    d = (A - K*G)*eps_hat + A*H*ukmin1 + K*shatk;
    uk = lsqlin(C,d);
    eps_hatkplus1 = (A - K*G)*eps_hat + A*H*ukmin1 - H * uk + K*shatk;
    shatk = G*eps_hat;
    eps_hat_mean_removed(:,i) = eps_hat - mean(eps_hat);
    sigma = var(eps_hat_mean_removed);
    eps_hat = eps_hatkplus1;
end
var_eps = mean(sigma);
