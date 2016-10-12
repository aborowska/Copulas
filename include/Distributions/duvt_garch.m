function dens = duvt_garch(x, mu, Sigma, df, GamMat)
%  location/scale univariate student t distribution 

%     c1 = gamma((df+1)/2); 
%     c2 = gamma(df/2);
    c0 = df+1;
    if c0 < 200
        c1 =  GamMat(floor(c0*50000));   
    else
        c1 = gamma(c0/2);
    end        
        
    if df < 200
        c2 = GamMat(round(df*50000)); % read the correponding value of th gamma
    else
        c2 = gamma(df/2);
    end
 
%     c3 = (pi*df)^(d/2);
%     c = c1/(c2*c3*sqrt(det(Sigma)));
%     c = c1/(c2*c3*tmp);
    c = c1/(c2*sqrt(pi*df*Sigma));
%     e = -(d+df)/2;
    e = -(df+1)/2;
%     dens = c*(1+(x-mu)*inv(Sigma)*(x-mu)'/df)^e;
    tmp = (x-mu)*(x-mu)/Sigma;
   
    dens = c*(1+tmp/df)^e;
end