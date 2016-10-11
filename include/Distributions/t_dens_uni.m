function d = t_dens_uni(x,nu)
    c1 = gamma((nu+1)/2);
    c2 = sqrt(nu*pi);
    c3 = gamma(nu/2);
    e = -(nu+1)/2;
    c = c1/(c2*c3);
    d = c.*(1 + (x.^2)./nu).^e;
end