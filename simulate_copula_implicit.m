clear all
close all

d = 2;
N = 1000; 

Rho = [1,0.5; 0.5,1];

df = 4;


%% Bivariate
% Normal
X_b_n = mvnrnd([0,0], Rho, N);

% Student's t
% t is normal variance mixture: X = sqrt(W)*Z, W indep Z, Z~N, nu/W~chi2(nu) [equiv.: W~IG(nu/2,nu/2)]
X_b_t4 = mvnrnd([0,0],Rho, N);
W = df./chi2rnd(df,N,1);
W = sqrt(W);
X_b_t4 = bsxfun(@times, X_b_t4, W);

% Elliptical implied copulas with correlation 0.5

figure(1)
subplot(3,2,1)
scatter(X_b_n(:,1),X_b_n(:,2),'.')
corrcoef(X_b_n) % 0.4921
subplot(3,2,2)
scatter(X_b_t4(:,1),X_b_t4(:,2),'.')
corrcoef(X_b_t4) % 0.5174

%% Margins
% Normal
X_n_marg = randn(N,d); 

% Student's t 
X_t4_marg = trnd(df,N,d);

%% Copulas: elliptical
% Normal
X_n = norminv(normcdf(X_n_marg)); 

% Student's t 
X_t4 = tinv(tcdf(X_t4_marg,df),df); 

figure(1)
subplot(3,2,3)
scatter(X_n(:,1),X_n(:,2),'.')
corrcoef(X_n) % 0.4921
subplot(3,2,4)
scatter(X_t4(:,1),X_t4(:,2),'.')
corrcoef(X_t4) % 0.5174

%% Copulas: non-elliptical
% Normal
X_meta_n = norminv(tcdf(X_t4_marg,df)); 

% Student's t 
X_meta_t4 = tinv(normcdf(X_n_marg),df); 

figure(1)
subplot(2,2,5)
scatter(X_meta_n(:,1),X_meta_n(:,2),'r.')
corrcoef(X_meta_n(all(isfinite(X_meta_n),2),:))

subplot(2,2,6)
scatter(X_meta_t4(:,1),X_meta_t4(:,2),'.')
corrcoef(X_meta_t4) 

