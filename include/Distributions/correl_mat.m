function Rho = correl_mat(rhos)
% create a dxd correlation matrix from a (d*(d-1)/2)x1 vector 
    [m,n] = size(rhos);
    rhos = reshape(rhos,m*n,1);
    m = m*n;
    d = (1 + sqrt(1+8*m))/2;
    if mod(d,1)~=0 % not an integer
        error('Input vecor of inccorect size.')
    else
        Rho = triu(ones(d),1)';
        Rho(Rho==1) = rhos;
        Rho = Rho + Rho' + diag(ones(d,1));
    end
end