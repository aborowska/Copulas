function  ep = duvt(eps, nu, hp, L)  
% density of the univariate t distribution
% calculated at matrix eps, with hp columns
% df of each row are in the nu vector
    c = gamma((nu+1)/2)./(sqrt(pi*nu).*gamma(nu/2));
    c = repmat(c, 1, hp);
    e = -(nu+1)/2;
    e = repmat(e, 1, hp);
    ep = c.*((1 + (eps.^2)./repmat(nu, 1, hp)).^e);

    if L
        ep = log(ep);
        ep = sum(ep,2);
    else
        ep = prod(ep,2);
    end
end