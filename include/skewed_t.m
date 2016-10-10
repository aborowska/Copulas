function val = skewed_t(z,nu,lambda)
% z = -2:0.01:2;
% nu = 3;
% lambda = 0;
C = gamma((nu+1)/2)/(gamma(nu/2)*sqrt(pi*(nu-2)));
A = 4*lambda*C*(nu-2)/(nu-1);
B = sqrt(1+3*lambda^2-A^2);

val = ((B.*z+A)./(1+sign(z+A/B).*lambda)).^2;
val = val./(nu-2);
val = B*C*(1+val).^(-(nu+1)/2);
end

figure(11)
plot(z,val)

hold on
plot(z,tpdf(z,nu),'r')
hold off
