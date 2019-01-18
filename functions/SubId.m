function [As, Cs, Ks] = SubId(s_id, N_id, N_val, s, n)

% Construction block hankel matrices Y0sN and YssN from s_id
Y0sN = zeros(72*s, (N_id));
for i=1:N_id
    Y0sN(:,i) = [s_id(:,i); s_id(:,i+1); s_id(:,i+2); s_id(:,i+3); s_id(:,i+4); s_id(:,i+5);
                s_id(:,i+6); s_id(:,i+7); s_id(:,i+8); s_id(:,i+s-1)];
end
YssN = Y0sN(:, s:end);

% Taking QR factorization of [Y0sN; YssN]
r = triu( qr( [Y0sN, YssN]' ) )';
R = r(:,1:720); 

% Taking SVD of ans
[U,S,V] = svd(inv(R)*YssN);

% Finding estimate of Xhat_sN
XhatsN = S.^(0.5) * V;

% Solving LSP
d = -[ XhatsN(:,s+1:end); YssN(:,1:end-1)];
C = XhatsN(:,s:end-1);
% minimizes x for norm(Cx-d), where x = [As, Cs]'
[x,resnorm,residual,exitflag,output,lambda]=lsqlin(C,d);

% Getting residuals:
As_res=residual(1:72,:);
Cs_res=residual(73:end,:);

% Calculating Covariance matrices of residuals
Qhat =;
Shat =;
Rhat =;

% Calculating Phat from dare 
Phat=dare(As,Cs,Qhat,Rhat,Shat);

% Calculating Ks
Ks = (Shat + As*Phat*Cs')(Rhat + Cs*phat*Cs')^(-1);