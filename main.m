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
phi_est     = pinv(G)*sk;

sumNom = 0;
sumDenom = 0;
phik = cell2mat(phiSim(1));
for i = 1:5000
    phik_mean_removed = phik(:,i) - mean(phik(:,i));
    phi_est_mean_removed = phi_est(:,i) - mean(phi_est(:,i));
    sumNom = sumNom + norm(phik_mean_removed - phi_est_mean_removed);
    sumDenom = sumDenom + norm(phik_mean_removed);
end
nom = sumNom/5000;
denom = sumDenom/5000;
VAFnocontrol = max(0,(1-nom/denom)*100);

%% 1.2 Unbiased Minimum Variance Estimate

usedPhiIdent = cell2mat(phiIdent(1));

% Calculating Covariance Matrix
summedPhi = 0;
for i = 1:length(usedPhiIdent)
    phiIdentMeanRemoved = usedPhiIdent(:,i) - mean(usedPhiIdent(:,i));
    summedPhi = summedPhi + phiIdentMeanRemoved*phiIdentMeanRemoved';
end
C_phi0     = (1/length(usedPhiIdent))*summedPhi;
% Expressing Unbiased Minimum Variance Estimate
phi_hat     = (inv(C_phi0) + G' * G)^(-1) * G' * sk;

%% 1.3
usedPhiSim = cell2mat(phiSim(1));
N = length(usedPhiSim);
residual_wave_front = zeros(49,N);
residual_wave_front(:,1) = usedPhiSim(:,1);
for i = 2:N
    residual_wave_front(:,i) = usedPhiSim(:,i) - phi_est(i-1);
end

%% 1.6 Random Walk function

% Calling function
[var_eps] = AOloopRW(G,H,C_phi0,sigmae,phiSim(1));

% Evaluating results
% if sigmaNoControl < var_eps
%     disp("it is better to not use the random walk method")
% elseif sigmaNoControl == var_eps
%     disp("using the random walk method yield neither better nor worse results")
% else
%     disp("using the random walk method is better than using no control")
% end
%% calculate VAF
sumNom = 0;
sumDenom = 0;
nom = 0;
denom = 0;
phik = cell2mat(phiSim(1));
for i = 1:5000
    phik_mean_removed = phik(:,i) - mean(phik(:,i));
    sumNom = sumNom + norm(phik_mean_removed + var_eps);
    sumDenom = sumDenom + norm(phik_mean_removed);
end
nom = sumNom/5000;
denom = sumDenom/5000;
VAFrw = max(0,(1-nom/denom)*100);
%% compute covariances

summedPhi = 0;
meanlessPhiSim = usedPhiSim - mean(usedPhiSim);
for i = 1:length(usedPhiSim)
    summedPhi = summedPhi + meanlessPhiSim(:,i)*meanlessPhiSim(:,i)';
end
C_phi0 = summedPhi/length(usedPhiSim);
summedPhi = 0;
for i = 2:length(usedPhiSim)
    summedPhi = summedPhi + meanlessPhiSim(:,i)*meanlessPhiSim(:,i-1)';
end
C_phi1 = summedPhi/(length(usedPhiSim)-1);
%% 2.7 compute kalman matrix and var_eps of VAR model

[ A, Cw, K ] = computeKalmanAR(C_phi0, C_phi1, G, sigmaNoControl);

[var_eps] = AOloopAR(G,H,C_phi0,sigmae,A,Cw,K,phiSim(1));


%% 3.5 Subspace Identification


s=10;           
n=500;           %order of system



% generating noisy wavefront signal s_id
phi_id = cell2mat(phiIdent(1));
s_id = G*phi_id + randn(72,5000)*sigmae^2; 

%% TEST QR
Y0sN = zeros(72*s, (n));
for i=1:n
    Y0sN(:,i) = [s_id(:,i); s_id(:,i+1); s_id(:,i+2); s_id(:,i+3); s_id(:,i+4); s_id(:,i+5);
                s_id(:,i+6); s_id(:,i+7); s_id(:,i+8); s_id(:,s)];
end
YssN = Y0sN(:, s:end);

r = triu(qr([Y0sN ; YssN]'))';
 

