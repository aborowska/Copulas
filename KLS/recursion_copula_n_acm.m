function [f,rho]=recursion_gam(y,PARAM,horizon)
%% Gaussian Copula GAS forecasting

n=size(y,2);

omega=PARAM(1);
A=PARAM(2);
B=PARAM(3);

f=zeros(1,n);
rho=zeros(1,n);
f(:,1)=(omega/(1-A-B));
rho(:,1)=f(:,1);
 for i=1:n-1
     s=y(1,i)*y(2,i);
     f(1,i+1)=omega+A*s+B*f(1,i);  
     if f(1,i+1)>=1
        f(1,i+1)=0.99999999;
     elseif f(1,i+1)<=-1
        f(1,i+1)=-0.9999999;
     end
     rho(1,i+1)=f(1,i+1);    
 end
 
 
f=f(1,end-horizon+1:end);
rho=rho(1,end-horizon+1:end);
 
end



