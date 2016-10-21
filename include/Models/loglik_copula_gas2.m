function [LL, f, rho] = loglik_copula_gas2(theta, y)
    [N,~] = size(theta);
    T = size(y, 1);
    
    omega = theta(:,1);
    A = theta(:,2);
    B = theta(:,3);
 
    z = norminv(y);
    
    f = zeros(N,T);  % alpha in the paper, time-varying parameter following GAS recursion
    f(:,1) = omega./(1-B); 
%     transf = @(aa) (1 - exp(-aa))./(1 + exp(-aa));
    rho = zeros(N,T); % the time-varying parameter tranfsormed by the trnasf function --> rho, i.e. the correlation
    rho(:,1) = transf(f(:,1));
    
%     scoref = @(z1, z2, r) ((1 + r.^2).*(z1.*z2 - r) - r.*(z1.^2 + z2.^2 - 2))./((1 - r.^2).^2);
%     inff = @(r) (1 + r.^2)./((1 - r.^2).^2);

    for jj = 2:T
        s = scoref(z(jj-1,1), z(jj-1,2), rho(:,jj-1));
        scscore = s./sqrt(inff(rho(:,jj-1)));
        f(:,jj) = omega + A.*scscore + B.*f(:,jj-1);
        rho(:,jj) = transf(f(:,jj));        
    end
    
    Z1 = repmat(z(:,1)',N,1);
    Z2 = repmat(z(:,2)',N,1);    
    LL = - 0.5*log(1-rho.^2) - 0.5*(Z1.^2 + Z2.^2 - 2*rho.*Z1.*Z2)./(1-rho.^2) ...
         + 0.5*Z1.^2 + 0.5*Z2.^2;
    LL = -sum(LL,2)/T;

%     f_T = f(:,T);
%     rho_T = rho(:,T);
end

function s = scoref(z1, z2, r) 
    s = ((1 + r.^2).*(z1.*z2 - r) - r.*(z1.^2 + z2.^2 - 2))./((1 - r.^2).^2);
end

function invf = inff(r) 
    invf = (1 + r.^2)./((1 - r.^2).^2);
end

function aa = transf(aa)
    aa = (1 - exp(-aa))./(1 + exp(-aa));
end