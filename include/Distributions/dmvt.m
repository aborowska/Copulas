function dens = dmvt(x, mu, Sigma, df, GamMat)
%  location/scale multivariate student t distribution 

    [~,d] = size(mu);
% % % % 	c1 = gamma((df+d)/2);
% % % %     c2 = gamma(df/2);
%%
%     global GamMat

%     tmp = (df+d)/2;
%     tmp = round((df+d)*100000);   % round to the 4th decimal number
%     c1 = GamMat(tmp); % read the correponding value of th gamma 
    c0 = df+d;
    if c0 < 100
        c1 =  GamMat(floor(c0*50000));   
    else
        c1 = gamma(c0/2);
    end
%     tmp = df/2;
%     tmp = round(tmp*100000);   
%     c2 = GamMat(tmp); % read the correponding value of th gamma 
     
    if df < 100
        c2 = GamMat(round(df*50000)); % read the correponding value of th gamma
    else
        c2 = gamma(df/2);
    end
        
%%
%     c3 = (pi*df)^(d/2);
%      c = c1/(c2*c3*sqrt(det(Sigma)));
    if d == 1
        tmp = sqrt(Sigma);
    else
        tmp = sqrt(det(Sigma));
    end
%     c = c1/(c2*c3*tmp);
    c = c1/(c2*((pi*df)^(d/2))*tmp);
%     e = -(d+df)/2;
    e = -c0/2;
%     dens = c*(1+(x-mu)*inv(Sigma)*(x-mu)'/df)^e;
    tmp = (x-mu)/Sigma;
    tmp = sum(tmp.*(x-mu));
    dens = exp(log(c) + log((1+tmp/df)^e));
end