function [theta, lnk, ind_red, x, lng_y, lnw_x,  eps_bar, eps_sim, C_T, lnp_T, RND] = fn_rmvgt_robust(N, mit, kernel, resampl_on, extended)
% robust sampling from mixture of multivariate t densities
% when resampl_on == 1 samples are redrawn from mit if they correspond to a bad region with zero
% kernel density  (i.e. these with -Inf weights)     
% "Standard" kernel evaluation is modified to account for the state space models 
% (in that case there are five inputs and additional outputs)
        
    % Mixtue of t (mit) parameters:
    mu = mit.mu;
    Sigma = mit.Sigma;
    df = mit.df;
    p = mit.p;

    theta = rmvgt2(N, mu, Sigma, df, p); % Sampling from the mixture of t

    % lnk - N vector of log-kernel evaluations at draws
    fprintf('\nKernel computation \n')

    if (nargin == 4) % not Extended version (an observation driven model)
        lnk = kernel(theta);
    else % Extedned version (a parametric model)
        [lnk, x, lng_y, lnw_x, eps_bar, eps_sim, C_T, lnp_T, RND] =  kernel(theta);
    end

    ind_red = 0;
    if resampl_on % Perform resampling to consturct the required mixture      
        while any(lnk == -Inf)
            ind_red = ind_red + 1;
            ind = find(lnk == -Inf);
            n_resamp = length(ind);
%             fprintf('resampling %d draws.\n', n_resamp)
            draw_new = rmvgt2(n_resamp, mu, Sigma, df, p);
            theta(ind,:) = draw_new;

            if (nargin == 4)  % not Extedned version (an observation driven model)
                lnk(ind) = kernel(draw_new);
            else % Extedned version (a parametric model)
%                 [lnk(ind), x(ind,:), lng_y(ind), lnw_x(ind), x_smooth(:,ind)] = kernel(draw_new);         
                [lnk(ind), x(ind,:), lng_y(ind), lnw_x(ind)] = kernel(draw_new);         
            end
        end
    end
end 