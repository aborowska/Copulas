function [theta, t] = rmvgt(N,mit)
% Random sampling from mixture of t densitites
% N - number of draws
% mit - list with parameters of mixture of t density
    tic
    [H,d] = size(mit.mu); % number of components, dimension of t distribution
    % sample membership
    memb = randsample(1:H,N,true,mit.p);
    % randsample(1:3,10,true,[0.1 0.3 0.6])
    theta = zeros(N, d);
    for h=1:H
        ind_h = (memb == h);
        n_h = sum(ind_h);
        if (n_h>0)
            mu_h = mit.mu(h,:);
            Sigma_h = mit.Sigma(h,:);
            Sigma_h = reshape(Sigma_h,d,d);
            df_h = mit.df(h);
%             draw_h = mvtrnd(Sigma_h,df_h,n_h);
            draw_h = rmvt(mu_h,Sigma_h,df_h,n_h);
%             draw_h = draw_h + repmat(mu_h,n_h,1);
            theta(ind_h,:) = draw_h;
        end
    end
    if (nargout > 1)
        t = toc;
    end
end
