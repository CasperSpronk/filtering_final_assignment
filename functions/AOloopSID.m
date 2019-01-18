function [var_eps] = AOloopSID(G, H, As, Cs, Ks, sigmae, phi_id)

% Taking Singular Value Decomposition of G
[U,S,V]     = svd(G,'econ');                                                   
% Partitioning SVD as in (Verhaegen, p.34)
r           = rank(G);
U1          = U(:,1:r);
S1          = S(1:r,1:r);
V1          = V(:,1:r);

% creating x(k)
x = XhatsN(1,:);

% finding u(k) using lsqlin
C = -H;
d = -(V1*S1*U1)*Cs*As*x;
[u,resnorm,residual,exitflag,output,lambda]=lsqlin(C,d);

% calculating variance of residual

    
end