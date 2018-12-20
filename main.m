addpath("./standard_files");
addpath("./functions");
[sigmaNoControl, sk] = AOloop_nocontrol(phiSim(1,4),sigmae,H,G);

%% 1.1 reconstructing wavefront
clc
[U,S,V] = svd(G,'econ');
r = rank(G);
S1 = S(1:r,1:r);
U1 = U(:,1:r);
V1 = V(:,1:r);
phi_est = V1*inv(S1)*U1'*sk;
F = V1*inv(S1)*U1';
((F'*F)^-1*F')'*sk;
F'*F;
r2 = rank(F'*F);