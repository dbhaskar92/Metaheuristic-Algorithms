% Generate random orthogonal matrix

orth(randn(5,5))

% Generate linear transformation matrix
% with condition number cn

nr = 5;
nc = 5;
cn = 3;
A=randn(nr,nc);
[U,S,V]=svd(A);
S(S~=0)=linspace(cn,1,min(nr,nc));
A=U*S*V'
cond(A)