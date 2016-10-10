function [PARAM,f,theta]=estimate_gam(y,PARAM)
% Gaussian Copula GAM estimation

%% Estimation
options=optimset('display','off','TolFun',1e-5,'LargeScale','off','TolX',1e-5,'HessUpdate','bfgs','FinDiffType','central',...
     'maxiter',5000,'MaxFunEvals',8000);
aux=PARAM;
PARAM(3)=invlogsig(PARAM(3));
PARAM(2)=invlogsig(PARAM(2)/(1-aux(3)));
PARAM=fminunc(@(PARAM) Lik(PARAM,y),PARAM,options);

%% Finalizing
[~,f,theta]=Lik(PARAM,y);
PARAM(3)=logsig(PARAM(3));
PARAM(2)=logsig(PARAM(2))*(1-PARAM(3));

end

function [L,f,rho]=Lik(PARAM,y)
%% Initialization
n=size(y,2);

omega=PARAM(1);
B=logsig(PARAM(3));
A=logsig(PARAM(2))*(1-B);


f=zeros(1,n);
f(:,1)=(omega/(1-A-B));
 for i=1:n-1
     if f(1,i)>=1
        f(1,i)=0.99999999;
     elseif f(1,i)<=-1
        f(1,i)=-0.9999999;
     end
     s=y(1,i)*y(2,i);
     f(1,i+1)=omega+A*s+B*f(1,i);       
end
rho=f;    

  

%% Likelihood
X=(y(1,:).^2+y(2,:).^2-2*rho.*y(1,:).*y(2,:))./(1-rho.^2);
L=-(1/2)*log(1-rho.^2)-1/2*X;
L=-sum(L,2)/n;

if isnan(L)==1||isinf(L)==1
   error('Likelihood error: GAM');
end

end


function y=invlogsig(x)
y=-log((1-x)./x);
end

