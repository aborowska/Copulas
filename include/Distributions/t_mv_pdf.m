function pdf = t_mv_pdf(x,mu,Sigma,nu)
    [~, d] = size(x);

    pdf = log(gamma((nu+d)/2)) - log(gamma(nu/2)) - 0.5*d*log(pi*nu);
    if isempty(mu)
        mu = zeros(1,d);
    end
    
    if (d == 1)
        pdf = pdf - 0.5*log(Sigma);
    elseif (d ==2)
        pdf = pdf - 0.5*log(Sigma(1,1)*Sigma(2,2) - Sigma(2,1)*Sigma(1,2));
    else
        pdf = pdf - 0.5*log(det(Sigma));
    end
    pdf = pdf - ((nu+d)/2)*log(1 + sum((bsxfun(@minus,x,mu)/Sigma).*bsxfun(@minus,x,mu),2)/nu);

    pdf = exp(pdf);
end