function [f,rho]=recursion_gas1(y,PARAM,horizon)
%% Student Copula GAS forecasting

n=size(y,2);

omega=PARAM(1);
A=PARAM(2);
B=PARAM(3);
nu=PARAM(4);

y=tinv(y,nu);

f=zeros(1,n);
rho=zeros(1,n);
f(:,1)=(omega/(1-B));
rho(:,1)=transf(f(:,1));
 for i=1:n-1
     s=score(y(:,i),rho(1,i),nu);
     f(1,i+1)=omega+A*s+B*f(1,i);  
     if f(1,i+1)>=5
        f(1,i+1)=5;
     elseif f(1,i+1)<=-5
        f(1,i+1)=-5;
     end
     rho(1,i+1)=transf(f(1,i+1));      
 end
 
 
f=f(1,end-horizon+1:end);
rho=rho(1,end-horizon+1:end);
 
end

function s=score(y,theta,nu)
%m=(y(1,1)^2+y(2,1)^2-2)/((1-theta^2)^2);
m=(y(1,1)^2+y(2,1)^2-2*theta*y(1,1)*y(2,1))/(1-theta^2);
omega=(nu+2)/(nu+m);
score=((1+theta^2)*(omega*y(1,1)*y(2,1)-theta)-theta*(omega*y(1,1)^2+omega*y(2,1)^2-2))/((1-theta^2)^2);
I=((1+theta^2-2*(theta^2)/(nu+2))/((1-theta^2)^2))*((nu+2)/(nu+4));
s=score/sqrt(I);
%s=score/I;
end

function theta=transf(f)
% for some reason using this instead of calling logsig()
% considerably speeds things up
theta=2*(1./(1+exp(-f))-0.5);
end
