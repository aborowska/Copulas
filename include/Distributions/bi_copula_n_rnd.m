function u_sim = bi_copula_n_rnd(rho)
% drawing from bivariate Gaussian copula
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
    
%     u_sim = normcdf(draw); 
    u_sim = 0.5 * erfc(-draw ./ sqrt(2));
end
