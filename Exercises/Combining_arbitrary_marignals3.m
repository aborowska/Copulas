% START WITH BIVARIATE T
clear all
close all

m = 2;
N = 2000;
Sigma = [1, 0.8; 0.8, 1];
C = chol(Sigma);
df1 = 3;
df2 = 10;

X = randn(N,m);
X = X*C;

W1 = df1./chi2rnd(df1,N,1);
X1 = bsxfun(@times,X,sqrt(W1));
U1 = tcdf(X1,df1); % t copula with df1 df

W2 = df2./chi2rnd(df2,N,1);
X2 = bsxfun(@times,X,sqrt(W1));
U2 = tcdf(X2,df2); % t copoula with df2 df

Y1 = zeros(N,m);
Y2 = zeros(N,m);
Y3 = zeros(N,m);
Y4 = zeros(N,m);
Y5 = zeros(N,m);
Y6 = zeros(N,m);

Y1(:,1) = tinv(U1(:,1),df1);
Y1(:,2) = tinv(U1(:,2),df1);

Y2(:,1) = tinv(U1(:,1),df1);
Y2(:,2) = tinv(U1(:,2),df2);

Y3(:,1) = tinv(U1(:,1),df2);
Y3(:,2) = tinv(U1(:,2),df2);

Y4(:,1) = tinv(U2(:,1),df1);
Y4(:,2) = tinv(U2(:,2),df1);

Y5(:,1) = tinv(U2(:,1),df1);
Y5(:,2) = tinv(U2(:,2),df2);

Y6(:,1) = tinv(U2(:,1),df2);
Y6(:,2) = tinv(U2(:,2),df2);

subplot(2,3,1)
scatter(Y1(:,1),Y1(:,2))
subplot(2,3,2)
scatter(Y2(:,1),Y2(:,2))
subplot(2,3,3)
scatter(Y3(:,1),Y3(:,2))
subplot(2,3,4)
scatter(Y4(:,1),Y4(:,2))
subplot(2,3,5)
scatter(Y5(:,1),Y5(:,2))
subplot(2,3,6)
scatter(Y6(:,1),Y6(:,2))