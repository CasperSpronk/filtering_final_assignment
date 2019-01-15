%% Delft University of Technology
%  SC42025: Filtering and Identification
%  Final Assignment; Turbulence Modeling for Adaptive Optics
%  By: Jorge   Bonekamp, 4474554
%      Casper  Spronk,   4369475
%  Date: 
%---------------------------------------------------------
clc
clear all


addpath("./standard_files");
addpath("./functions");
load('systemMatrices.mat');
load('turbulenceData.mat');

[sigmaNoControl, sk] = AOloop_nocontrol(phiSim(1),sigmae,H,G);

%% 1.1 Reconstructing wavefront
clc
% Taking Singular Value Decomposition of G
[U,S,V]     = svd(G,'econ');                                                   
% Partitioning SVD as in (Verhaegen, p.34)
r           = rank(G);
U1          = U(:,1:r);
S1          = S(1:r,1:r);
V1          = V(:,1:r);
% Expressing estimate of the wavefront
phi_est     = V1*inv(S1)*U1'*sk;

% F = V1*inv(S1)*U1';
% ((F'*F)^-1*F')'*sk;
% F'*F;
% r2 = rank(F'*F);

%% 1.2 Unbiased Minimum Variance Estimate

usedPhiIdent = cell2mat(phiIdent(1));
% Calculating Covariance Matrix
summedPhi = 0;
for i = 1:length(usedPhiIdent)
    summedPhi = summedPhi + usedPhiIdent(:,i)*usedPhiIdent(:,i)';
end
C_phi     = (1/length(usedPhiIdent))*summedPhi;
% Expressing Unbiased Minimum Variance Estimate
phi_hat     = (inv(C_phi) + G' * G)^(-1) * G' * sk;

%% 1.3
usedPhiSim = cell2mat(phiSim(1));
N = length(usedPhiSim);
residual_wave_front = zeros(49,N);
residual_wave_front(:,1) = usedPhiSim(:,1);
for i = 2:N
    residual_wave_front(:,i) = usedPhiSim(:,i) - phi_est(i-1);
end

%% 1.6 Random Walk function
[var_eps] = AOloopRW(G,H,C_phi,sigmae,phiSim(1));

if sigmaNoControl < var_eps
    disp("it is better to not use the random walk method")
elseif sigmaNoControl == var_eps
    disp("using the random walk method yield neither better nor worse results")
else
    disp("using the random walk method is better than using no control")
end