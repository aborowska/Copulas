function [PARAM,f,theta]=estimate_gas1(y,PARAM)
% Gaussian Copula GAS estimation
% Time varying variance GAS estimation
% GAS 1: with parameter tranformation, square-root scaling
% GAS 2: with parameter tranformation, scaling by the inverse of the information matrix
% GAS 3: no parameter tranformation, square-root scaling
% GAS 4: no parameter tranformation, scaling by the inverse of the
% information matrix


%% Estimation
options=optimset('display','off','TolFun',1e-5,'LargeScale','off','TolX',1e-5,'HessUpdate','bfgs','FinDiffType','central',...
     'maxiter',10000,'MaxFunEvals',20000);
PARAM(3)=invlogsig(PARAM(3));
try
PARAM=fminunc(@(PARAM) Lik(PARAM,y),PARAM,options);
catch
PARAM=fminsearch(@(PARAM) Lik(PARAM,y),PARAM,options);
end

%% Finalizing
[~,f,theta]=Lik(PARAM,y);
PARAM(3)=logsig(PARAM(3));

end

function [L,f,rho]=Lik(PARAM,y)
%% Initialization
n=size(y,2);

omega=PARAM(1);
A=PARAM(2);
B=logsig(PARAM(3));

f=zeros(1,n);
rho=zeros(1,n);
f(:,1)=(omega/(1-B));
rho(:,1)=transf(f(:,1));
 for i=1:n-1
     s=score(y(:,i),rho(1,i));
     f(1,i+1)=omega+A*s+B*f(1,i);  
     if f(1,i+1)>=5
        f(1,i+1)=5;
     elseif f(1,i+1)<=-5
        f(1,i+1)=-5;
     end
     rho(1,i+1)=transf(f(1,i+1));    
 end


%% Likelihood
X=(y(1,:).^2+y(2,:).^2-2*rho.*y(1,:).*y(2,:))./(1-rho.^2);
L=-(1/2)*log(1-rho.^2)-1/2*X;
L=-sum(L,2)/n;

if isnan(L)==1||isinf(L)==1
   error('Likelihood error: GAS1');
end

end

function s=score(y,theta)
score=((1+theta^2)*(y(1,1)*y(2,1)-theta)-theta*(y(1,1)^2+y(2,1)^2-2))/((1-theta^2)^2);
I=(1+theta^2)/((1-theta^2)^2);
s=score/sqrt(I);
end

function theta=transf(f)
% for some reason using this instead of calling logsig()
% considerably speeds things up
theta=2*(1./(1+exp(-f))-0.5);
end

function y=invlogsig(x)
y=-log((1-x)./x);
end