function u_sim = bi_copula_t_rnd(rho, df)
% drawing from bivariate Student's copula
% with time varying correlation rho
    T = size(rho,2);
    Rho = ones(2,2,T);
    Rho(1,2,:) = rho;
    Rho(2,1,:) = rho;
    draw = mvnrnd([0,0],Rho);
    % SIGMA is a  D-by-D symmetric positive semi-definite matrix, 
    % or a D-by-D-by-N array.
    % If SIGMA is an array, MVNRND generates each row of draw using the
    % corresponding page of SIGMA, i.e., MVNRND computes R(I,:) using MU(I,:)
    % and SIGMA(:,:,I).
    R = df./chi2rnd(df,T,1);
    R = sqrt(R);
    draw = bsxfun(@times, draw, R);
    
    u_sim = tcdf(draw,df); 
end
