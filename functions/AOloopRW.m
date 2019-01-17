function [var_eps] = AOloopRW(G,H,C_phi_zero,sigmae,phi_sim)

[ sigma, sko ] = AOloop_nocontrol(phi_sim,sigmae,H,G);

usedPhiSim = cell2mat(phi_sim);
[n, m] = size(usedPhiSim);

[U,S,V]     = svd(G,'econ');
r           = rank(G);
U1          = U(:,1:r);
S1          = S(1:r,1:r);
V1          = V(:,1:r);

% reconstructing wavefront
phik_hat = V1*inv(S1)*U1'*sko;

phiKPlus1 = zeros(n,m);
uk = zeros(n,m);
epsilon = zeros(n,m+1);
[mG, nG] = size(G);
sk = zeros(mG,m);
sigma = zeros(m,1);
epsilonkplus1k = zeros(n,m);
% disp(wgn(1,n,0))
% calculating sk

for i = 2:m
    phiKPlus1(:,i) = phik_hat(:,i-1) + randn(n,1);
%     uk(:,i) = H^-1 * phiKPlus1(:,i);
%     epsilon(:,i+1) = phiKPlus1(:,i) - H * uk(:,i);
%     sk(:,i) = G * epsilon(:,i) + randn(mG,1);
end
C = -H;
for i = 1:m-1
    if i == 1
        epsilonkplus1k(:,i+1) = (C_phi_zero^-1 + G'*G)^-1 * G'*sko(:,1);
        d = -(C_phi_zero^-1 + G'*G)^-1 * G'*sko(:,1);
    else
        epsilonkplus1k(:,i+1) = (C_phi_zero^-1 + G'*G)^-1 * G'*skplus1 - H*deltaUK;
        epsilonkplus1k(:,i+1) = (C_phi_zero^-1 + G'*G)^-1 * G'*skplus1 - H*deltaUK - mean(phik_hat(:,i));
        d = -(C_phi_zero^-1 + G'*G)^-1 * G'*skplus1;
    end
    deltaUK = lsqlin(C,d);   
    skplus1 = G*(phik_hat(:,i) + rand(n,1) - H*deltaUK) + randn(mG,1);
    phik_hat(:,i+1) = V1*inv(S1)*U1'*skplus1;
    
end
var_eps = var(epsilonkplus1k);
% var_eps = mean(var_eps);

% 
% epsilonHatKK = zeros(n,m);
% % calculating epsilon hat KK
% for i = 2:m
%     epsilonHatKK(:,i) = (C_phi_zero^-1 + G'*G)^-1 * G' * sk(:,i);
% end
% % 
% % C = -H;
% % deltaUK = zeros(n,m);
% % solving linear least squares problem formulated in 1.5
% % for i = 1:m
% %     d = epsilonHatKK(:,i);
% %     deltaUK(:,i) = lsqlin(C,d);
% % end
% HdeltaUK = H*deltaUK;
% sigma = zeros(m,1);
% epsilonHatK1K = zeros(n,m);
% % calculate epsilon hat (k+1|k) and variance of epsilon hat (k+1|k)
% for i = 1:m-1
%     epsilonHatK1K(:,i+1) = epsilonHatKK(i) - H * deltaUK(:,i);
%     sigma(m) = var(epsilonHatK1K(:,i+1));
% end
% 
% var_eps = mean(sigma);



