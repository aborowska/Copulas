function loglik = loglik_t_copula(u,theta)
    N = size(u,1);
    nu = theta(:,1);
    rho = theta(:,2);
    
    % transform parameters    
    nu = 2 + exp(nu);
    rho = exp(rho)/(1+exp(rho));
    
    u = tinv(nu,nu);
    
    % t-copula
    t_copula_pdf = log(gamma((nu+2)/2)) + log(gamma(nu/2)) - 2*log(gamma((nu+1)/2))...
    - 0.5*log(1 - rho^2)...
    - ((nu+2)/2)*log(1 + (u(:,1).^2 + u(:,2).^2 - 2*rho*u(:,1).*u(:,2))/(nu*(1-rho^2))) ...
    + ((nu+1)/2)*log(1 + (u(:,1).^2)/nu) ...
    + ((nu+1)/2)*log(1 + (u(:,2).^2)/nu);   
   
    % average minus loglik
    loglik = -sum(t_copula_pdf)/N;

end