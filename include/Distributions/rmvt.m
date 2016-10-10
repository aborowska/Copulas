function draw = rmvt(mu,Sigma,df,n)
%  location/scale multivariate student t distribution 
    [N_mu,d] = size(mu);
%     mu_n = zeros(1,d);
%     Z = mvnrnd(mu_n,Sigma,n);
% %     R2 = chi2rnd(df,n,1);
% %     R2 = R2/df;
% %     R2 = repmat(R2,1,d);
% % %     draw = repmat(mu,n,1) + Z./sqrt(R2);
% %     R2 = chi2rnd(df,n,1);
% %     df = df*ones(n,1);
% %     R2 = df./R2;
%     R2 = random('gam',df/2,2/df,n,1); % shape scale
%     R2 = repmat(R2,1,d);
%     draw = repmat(mu,n,1) + Z.*sqrt(R2);

%     R = chol(Sigma);
%     Y = mvtrnd(eye(d), df, n);
    Y = mvnrnd(zeros(1,d),Sigma,n);
    R = df./chi2rnd(df,n,1);
    R = repmat(R,1,d);
    if (N_mu == 1)
        draw = repmat(mu,n,1) + Y.*sqrt(R);
    else
        draw = mu + Y.*sqrt(R);
    end
end