function [f,rho]=recursion_gas1(y,PARAM,horizon)
%% Gaussian Copula GAS forecasting

n=size(y,2);

omega=PARAM(1);
A=PARAM(2);
B=PARAM(3);


f=zeros(1,n);
rho=zeros(1,n);
f(:,1)=(omega/(1-B));
rho(:,1)=transf(f(:,1));
 for i=1:n-1
     s=score(y(:,i),rho(1,i));
     f(1,i+1)=omega+A*s+B*f(1,i);  
     rho(1,i+1)=transf(f(1,i+1));      
 end
 
f=f(1,end-horizon+1:end);
rho=rho(1,end-horizon+1:end);
 
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
