clear all
close all
N = 2000;
d = 2;
rho = 0.8;
Sigma = correl_mat(rho);
df = 5;
gam = [0.8 0];

C = chol(Sigma);
X = randn(N,d);
X = X*C;

W = df./chi2rnd(df,N,1);
X = bsxfun(@times,X,sqrt(W));
Xsk = bsxfun(@plus,X,gam*W);
hold on
scatter(X(:,1),X(:,2),'.')
scatter(Xsk(:,1),Xsk(:,2),'r.')
hold off