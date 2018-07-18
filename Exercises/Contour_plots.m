% contour plots
M = 100;
OO = ones(1,M);
xx = linspace(-4,4,M)';
vv = (0.02:0.02:0.28);

t_uv_pdf = @(x,mu,Sigma,nu) exp(log(gamma((nu+1)/2)) - log(gamma(nu/2)) - 0.5*log(pi*nu) - 0.5*log(Sigma)...
    - ((nu+1)/2)*log(1 + ((x-mu).^2)/(Sigma*nu)));
t_copula_pdf = @(u,R,nu) t_mv_pdf(tinv(u,nu),[],R,nu)./prod(t_uv_pdf(tinv(u,nu),0,1,nu),2);

rho = 0.8;
Rho = Correl_mat(rho);

DF = [3,10];
figure(1)
for pp = 1:6
    % margins:
    if (mod(pp,3) == 0)
        df1 = DF(2);
    else
        df1 = DF(1);
    end
    if (mod(pp,3) == 1)
        df2 = DF(1);
    else
        df2 = DF(2);
    end
    xx1 = tpdf(xx,df1);
    xx2 = tpdf(xx,df2);

    % feeds to copula:
    uu1 = tcdf(xx,df1);
    uu2 = tcdf(xx,df2);

    % copula common evaualtion
    if (pp <= 3)
        df = DF(1);
    else
        df = DF(2);
    end
    zz1 = tinv(uu1,df);
    zz2 = tinv(uu2,df);

    tt1 = t_uv_pdf(zz1,0,1,df);
    tt2 = t_uv_pdf(zz2,0,1,df);

    % copula pdf
    CPDF = zeros(M,M);
    for ii = 1:M
        u_in = [uu1(ii)*OO', uu2];
        CPDF(:,ii) = t_copula_pdf(u_in,Rho,df);
        % uu's are distinct feeds to copula
        % z_in is their common evaluation
    %     z_in = [tinv(uu1,df),tinv(uu2(ii),df)*ones(M,1)];
    %     CPDF(:,ii) = t_mv_pdf(z_in,[],Rho,df);
    %     CPDF(:,ii) = CPDF(:,ii)./prod(t_uv_pdf(z_in,0,1,df),2);
    % %     CPDF(:,ii) = CPDF(:,ii)./(tt1*tt2(ii));
    end
    % surf(CPDF)
    % contour(xx,xx,CPDF)

    % joint pdf = prod(margins pdf)*copula pdf
    Y = (OO'*xx1') .* (xx2*OO) .* CPDF;
    % surf(Y)
    subplot(2,3,pp)
    contour(xx,xx,Y,vv);
    title(sprintf('df1 = %i, df2 = %i, df = %i ',df1,df2,df));
end