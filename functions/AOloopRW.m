function [var_eps] = AOloopRW(G,H,C_phi_zero,sigmae,phi_sim)

[ sigma, sko ] = AOloop_nocontrol(phi_sim,sigmae,H,G);


phi_hat = (inv(C_phi_zero) + G' * G)^(-1) * G' * sk;

usedPhiSim = cell2mat(phi_sim);
% Taking Singular Value Decomposition of G
[U,S,V]     = svd(G,'econ');                                                   
% Partitioning SVD as in (Verhaegen, p.34)
r           = rank(G);
U1          = U(:,1:r);
S1          = S(1:r,1:r);
V1          = V(:,1:r);
% Expressing estimate of the wavefront
phi_est     = V1*inv(S1)*U1'*sko;

% calculating residual wave front
N = length(usedPhiSim);
residual_wave_front = zeros(49,N);
residual_wave_front(:,1) = usedPhiSim(:,1);
for i = 2:N
    residual_wave_front(:,i) = usedPhiSim(:,i) - phi_hat(i-1);
end

% calculate mean of residual wave front
residual_wave_front_mean = mean(residual_wave_front);
residual_wave_front_removed = residual_wave_front - residual_wave_front_mean;

n = size(H,1);      % dimension lifted wavefront
ns = size(G,1);     % dimension lifted sensor slopes
T = length(phik);   % number of temporal phase points

skRW = zeros(ns,T);   % slopes measurements
var_eps = zeros(T,1);
epsk = zeros(n,T);
eps_piston_removed = zeros(n,T); % residual wavefront with mean removed


for k = 2:T-1
    epsk(:,k) = phi_hat(:,k) - residual_wave_front_removed(:,k-1);
    eps_piston_removed(:,k+1) = epsk(:,k+1)-mean(epsk(:,k+1)); 
    skRW(:,k+1) = G*epsk(:,k+1) + sigmae*randn(ns,1);
    var_eps(k+1) = var(eps_piston_removed(:,k+1));
end
var_eps = mean(var_eps);





