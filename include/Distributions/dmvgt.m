function dens = dmvgt(theta, mit, L, GamMat)
% density of a mixture of multivariate t distributions
% L (log) - return log-density values if L=true
% % MATLAB
%     [H,d] = size(mit.mu); % number of components, dimension of t distribution
%     [N,~] = size(theta);
%     dcoms = zeros(N,H);
%     for h = 1:H
%         mu_h = mit.mu(h,:);
%         Sigma_h = mit.Sigma(h,:);
%         Sigma_h = reshape(Sigma_h,d,d);
%         df_h = mit.df(h);
%         for ii = 1:N
%             dcoms(ii,h) =  dmvt(theta(ii,:), mu_h, Sigma_h, df_h, GamMat);
%         end
%     end
%     tmp = log(repmat(mit.p,N,1)) + log(dcoms);
%     dens = sum(exp(tmp),2);
%     if (L == true)
%         dens = log(dens);
%     end
 
% % C-MEX
    L = double(L);
    dens = dmvgt_mex(theta, mit.mu, mit.Sigma, mit.df, mit.p, GamMat, L);

% % comparison of the MATLAB and C-MEX kernel evaluation functions
%     lnd_mex = dmvgt_mex(theta, mit.mu, mit.Sigma, mit.df, mit.p, GamMat, 1);
%     fprintf('\n *** sum(abs(lnd_mex-lnd)> eps)  = %6.4f ***\n', sum(abs(lnd_mex-dens)>eps) );
%     fprintf(' *** sum(abs(lnd_mex-lnd)) = %16.14f ***\n\n', sum(abs(lnd_mex-dens)));
%     if (sum(abs(lnd_mex-dens)) > 1e-4 )
%         keyboard
%         [val, MM] = max(abs(lnd_mex-dens));  % lnd_mex(MM)-dens(MM)
%     end
    
end