%% Delft University of Technology
%  SC42025: Filtering and Identification
%  Final Assignment; Turbulence Modeling for Adaptive Optics
%  By: Jorge   Bonekamp, 4474554
%      Casper  Spronk,   44
%  Date: 
%---------------------------------------------------------



addpath("./standard_files");
addpath("./functions");
load('systemMatrices.mat');
load('turbulenceData.mat');

%[sigmaNoControl, sk] = AOloop_nocontrol(phiSim(1,4),sigmae,H,G);

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


F = V1*inv(S1)*U1';
((F'*F)^-1*F')'*sk;
F'*F;
r2 = rank(F'*F);

%% 1.2 Unbiased Minimum Variance Estimate

% Calculating Covariance Matrix
% C_phi     = (1/N)*sum(  ;
% Expressing Unbiased Minimum Variance Estimate
phi_hat     = ( inv(C_phi) + G'G)^(-1)   G'sk;

