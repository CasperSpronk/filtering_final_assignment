function [As,Cs,Ks]=SubId(s_id, s, n)

% Construction block hankel matrices Y0sN and YssN
Y0sN = zeros(72*s, (n));
for i=1:n
    Y0sN(:,i) = [s_id(:,i); s_id(:,i+1); s_id(:,i+2); s_id(:,i+3); s_id(:,i+4); s_id(:,i+5);
                s_id(:,i+6); s_id(:,i+7); s_id(:,i+8); s_id(:,s)];
end
YssN = Y0sN(:, s:end);

% Taking QR factorization of [Y0sN; YssN]





% Taking SVD of ans




% Finding estimate of Xhat_sN






% Solving LSP