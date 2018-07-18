clear all
close all
addpath(genpath('include/'));

d = 4; % Dimension
partition = [1,3]; % subsets indicators
S = length(partition); % no of groups
df = [5,10]; % df of groups


N = 2000;
% rhos = [0.8 0.7 0.6 0.4 0.3 0.2];
rhos = [0.8 0.8 0.8 0.8 0.8 0.8];
Sigma = correl_mat(rhos);
C = chol(Sigma);

X = randn(N,d);
X = X*C; % sample from N(0,Sigma)
U = zeros(N,d);

Uw = rand(N,1);
W = zeros(N,S);

for s = 1:S
    W(:,s) = df(s)./chi2inv(Uw,df(s));
end 

for s = 1:S
    [s1, s2] = fn_partition_ends(partition, d, s);
    X(:,s1:s2) = bsxfun(@times,X(:,s1:s2),sqrt(W(1,s)));
    U(:,s1:s2) = tcdf(X(:,s1:s2),df(s));
end

Z = norminv(U);

subplot(4,4,1)
hist(Z(:,1))

subplot(4,4,6)
hist(Z(:,2))

subplot(4,4,11)
hist(Z(:,3))

subplot(4,4,16)
hist(Z(:,4))

%% Estimate t copula
    N = 2000;
    rho = 0.2;
    Rho = correl_mat(rho);
    C = chol(Rho);
    df = 3;

    X = randn(N,2);
    X = X*C;

    W = df./chi2rnd(df,N,1);
    X = bsxfun(@times,X,sqrt(W));
    U = tcdf(X,df); % t copula with df
    Z = norminv(U);
    u = norcdf(Z);

    theta = [10,0.5];
    theta(:,1) = log(theta(:,1) - 2);
    theta(:,2) = log(theta(:,2)/(1-theta(:,2)));
    
    f_opt = @(xx) loglik_t_copula(u,xx);
    f_opt(theta)
        
    theta = fminunc(f_opt,theta);
    
    nu = 2 + exp(theta(1,1));
    rho = exp(theta(1,2))/(1+exp(theta(1,2)));