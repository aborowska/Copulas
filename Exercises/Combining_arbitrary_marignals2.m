clear all
close all

m = 2;
N = 5000;
Sigma = [1, 0.5; 0.5, 1];
C = chol(Sigma);
df = 4;

% % Important: multiply a vector by the lower triangular!
% x = randn(m,1);
% x = C'*x; %!!!

%%
Z = randn(N,m);
Z = Z*C;

W = df./chi2rnd(df,N,1);
X = bsxfun(@times,Z,sqrt(W));


%% 
Z_pit = normcdf(Z); % normal copula
X_pit = tcdf(X,df); % t copula

corrcoef(Z)
corrcoef(X)
corrcoef(Z_pit)
corrcoef(X_pit)

scatter(Z(:,1),Z(:,2))
scatter(X(:,1),X(:,2))
scatter(Z_pit(:,1),Z_pit(:,2))
scatter(X_pit(:,1),X_pit(:,2))

%%
X_meta_n = tinv(Z_pit,df); % t marginas + normal copula = X_meta_n
Z_meta_t = norminv(X_pit); % n marginas + t copula = Z_meta_t

%%
figure(1)
subplot(2,2,1)
scatter(Z(:,1),Z(:,2),'.')
subplot(2,2,2)
scatter(X_meta_n(:,1),X_meta_n(:,2),'g.')
subplot(2,2,3)
scatter(Z_meta_t(:,1),Z_meta_t(:,2),'r.')
subplot(2,2,4)
scatter(X(:,1),X(:,2),'.')


figure(2)
subplot(2,3,1)
scatter(X_meta_n(:,1),X_meta_n(:,2),'g.')
subplot(2,3,2)
hist(X_meta_n(:,1))
subplot(2,3,3)
hist(X_meta_n(:,2))


subplot(2,3,4)
scatter(Z_meta_t(:,1),Z_meta_t(:,2),'r.')
subplot(2,3,5)
hist(Z_meta_t(:,1))
subplot(2,3,6)
hist(Z_meta_t(:,2))


