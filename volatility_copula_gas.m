function [f, rho] = volatility_copula_gas(theta, y)
    [N,~] = size(theta);
    T = size(y, 1);
    
    omega = theta(:,1);
    A = theta(:,2);
    B = theta(:,3);
 
    z = norminv(y);
    
    f = zeros(N,T);  % alpha in the paper, time-varying parameter following GAS recursion
    f(:,1) = omega./(1-B); 
    transf = @(aa) (1 - exp(-aa))./(1 + exp(-aa));
    rho = zeros(N,T); % the time-varying parameter tranfsormed by the trnasf function --> rho, i.e. the correlation
    rho(:,1) = transf(f(:,1));
    
    scoref = @(z1, z2, r) ((1 + r.^2).*(z1.*z2 - r) - r.*(z1.^2 + z2.^2 - 2))./((1 - r.^2).^2);
    inff = @(r) (1 + r.^2)./((1 - r.^2).^2);

    for jj = 2:T
        s = scoref(z(jj-1,1), z(jj-1,2), rho(:,jj-1));
        scscore = s./sqrt(inff(rho(:,jj-1)));
        f(:,jj) = omega + A.*scscore + B.*f(:,jj-1);
        rho(:,jj) = transf(f(:,jj));
    end
%     f_T = f(:,T);
%     rho_T = rho(:,T);
end

% function scaled_score = scaled_score_copula_gas(y1, y2, rho)
%     score = ((1 + rho.^2).*(y1.*y2 - rho) - rho.(y1.^2 + y2.^2 - 2))./((1 - rho.^2).^2);
%     I = (1 + rho.^2)./((1 - rho.^2).^2);
%     scaled_score = score./sqrt(I);
% end