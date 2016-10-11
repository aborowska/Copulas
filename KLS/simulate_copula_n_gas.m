function [y,rho]=simulate_gas(PARAM,n)
%% Gaussian Copula GAS simulate

omega=PARAM(1);
A=PARAM(2);
B=PARAM(3);

n=n+1000;
y=zeros(2,n);
f(:,1)=omega/(1-B);
f=repmat(f,1,n);
rho=repmat(transf(f(:,1)),1,n);
 for i=1:n-1
     U=copularnd('Gaussian',rho(1,i),1);
     y(1,i)=U(1,1);
     y(2,i)=U(1,2);
     y2=norminv(y(:,i));
     s=score(y2,rho(1,i));
     f(1,i+1)=omega+A*s+B*f(1,i);  
     rho(1,i+1)=transf(f(1,i+1));      
 end
U=copularnd('Gaussian',rho(1,end),1);
y(1,end)=U(1,1);
y(2,end)=U(1,2);
y=y(:,1001:end);
rho=rho(1,1001:end);
 
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
