%http://datascienceplus.com/modelling-dependence-with-copulas/
clear all
close all

m = 3;
N = 2000;
Sigma = [1, 0.4, 0.2; 0.4, 1, -0.8; 0.2, -0.8, 1];
C = chol(Sigma);

% Important: multiply a vector by the lower triangular!
x = randn(m,1);
x = C'*x; %!!!

%%
X = randn(N,m);
X = X*C;

RX = corrcoef(X);

figure(1)
subplot(3,3,1)
hist(X(:,1))
subplot(3,3,5)
hist(X(:,2))
subplot(3,3,9)
hist(X(:,3))
subplot(3,3,4)
scatter(X(:,1),X(:,2),'.')
subplot(3,3,7)
scatter(X(:,1),X(:,3),'.')
subplot(3,3,8)
scatter(X(:,2),X(:,3),'.')
subplot(3,3,2)
set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[]);
text(0.35,0.5,sprintf('%3.2f',RX(1,2)),'FontSize',14)
subplot(3,3,3)
set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[]);
text(0.35,0.5,sprintf('%3.2f',RX(1,3)),'FontSize',14)
subplot(3,3,6)
set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[]);
text(0.35,0.5,sprintf('%3.2f',RX(2,3)),'FontSize',14)

figure(6)
scatter3(X(:,1),X(:,2),X(:,3))

%% PITs
% And now comes the magic trick: recall that if X is a random variable with
% distribution F then F(X) is uniformly distributed in the interval [0, 1]. 
% In our toy example we already know the distribution F for each of the three random variables
% so this part is quick and straightforward.
U = normcdf(X);
RU = corrcoef(U);

figure(2)
subplot(3,3,1)
hist(U(:,1))
subplot(3,3,5)
hist(U(:,2))
subplot(3,3,9)
hist(U(:,3))
subplot(3,3,4)
scatter(U(:,1),U(:,2),'.')
subplot(3,3,7)
scatter(U(:,1),U(:,3),'.')
subplot(3,3,8)
scatter(U(:,2),U(:,3),'.')
subplot(3,3,2)
set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[]);
text(0.35,0.5,sprintf('%3.2f',RU(1,2)),'FontSize',14)
subplot(3,3,3)
set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[]);
text(0.35,0.5,sprintf('%3.2f',RU(1,3)),'FontSize',14)
subplot(3,3,6)
set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[]);
text(0.35,0.5,sprintf('%3.2f',RU(2,3)),'FontSize',14)

figure(3)
scatter3(U(:,1),U(:,2),U(:,3),'.')

%% Arbitrary marginals
% The last step: we only need to select the marginals and apply them to u.
% E.g. chose the marginals to be Gamma, Beta and Student distributed with the parameters specified below.

Y = U;
RY = corrcoef(Y);

Y(:,1) = gaminv(Y(:,1),2,1);
% X = gaminv(P,A,B) computes the inverse of the gamma cdf 
% with shape parameters in A and scale parameters in B 
Y(:,2) = betainv(Y(:,2),2,2);
Y(:,3) = tinv(Y(:,3),5);
figure(4)
scatter3(Y(:,1),Y(:,2),Y(:,3))

%% What is worth noticing is that by starting from a multivariate normal sample 
% we have build a sample with the desired and fixed dependence structure
% and, basically, arbitrary marginals.

figure(5)
subplot(3,3,1)
hist(Y(:,1))
subplot(3,3,5)
hist(Y(:,2))
subplot(3,3,9)
hist(Y(:,3))
subplot(3,3,4)
scatter(Y(:,1),Y(:,2),'.')
subplot(3,3,7)
scatter(Y(:,1),Y(:,3),'.')
subplot(3,3,8)
scatter(Y(:,2),Y(:,3),'.')
subplot(3,3,2)
set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[]);
text(0.35,0.5,sprintf('%3.2f',RY(1,2)),'FontSize',14)
subplot(3,3,3)
set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[]);
text(0.35,0.5,sprintf('%3.2f',RY(1,3)),'FontSize',14)
subplot(3,3,6)
set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[]);
text(0.35,0.5,sprintf('%3.2f',RY(2,3)),'FontSize',14)
