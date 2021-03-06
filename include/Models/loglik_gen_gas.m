function [LL, f, sigma2] = loglik_gen_gas(theta, y, link, scale, GamMat)
    % theta is Nxk, matrix of draws
    % k == 4 --> normal errors
    % k == 5 --> stundent's t errors, with nu degrees of freedom
    % (l)ink: 1 - linear link, sigma2(f)=f; 0 - exp link, sigma2(f)=exp(f); 
    % if exp link, nontrivial chain rule required, dsigma2/df = exp(f) = sigma2: 
    % (s)cale: 1 - inv fisher; 0 - sqrt inv fisher
    % Note: (l,s) == (1,1) was discussed in the original GAS paper (Creal et al. 2013)
    % Note: (l,s) == (0,0) cancels out so the updating step has a constant unit variance
    % and is invariant under any nondegenerate parameter transformation (cf. Koopman et al. 2016)
    [N, k] = size(theta);
    mu = theta(:,1);
    omega = theta(:,2);
    A = theta(:,3);
    B = theta(:,4); 
    
    if (k == 5)
        nu = theta(:,5);
%         rho = (nu-2)./nu;    
    %     nu_con1 = (nu+1)./(nu-2);
        nu_con2 = 2*((nu+3)./nu);
    end
    
    if link  
        fn_link = @(xx) xx;
        fn_chain_rule = @(xx) 1;
    else
        fn_link = @(xx) exp(xx);    
        fn_chain_rule = @(xx) xx;
    end
    
    if scale
        fn_scale = @(xx) xx;
    else
        fn_scale = @(xx) sqrt(xx);       
    end
     
    T = size(y,1);
     
%     d = -Inf*ones(N,1);
    LL = -Inf*ones(N,1);

    for ii = 1:N
        if mod(ii,1000) == 0
            fprintf('loglik ii = %d\n',ii);
        end
        
        f = zeros(T,1);
        sigma2 = zeros(T,1);
        
%         pdf = zeros(T,1);
        L = zeros(T,1);

        f(1,1) = omega(ii,1)/(1-B(ii,1)); % unconditional variance to initialize f_1
        sigma2(1,1) = fn_link(f(1,1));

%         pdf(1,1) = duvt_garch(y(1,1), mu(ii,1), rho(ii,1)*sigma2(1,1), nu(ii,1), GamMat);
%         pdf(1,1) = log(pdf(1,1));
%         y2 = ((y(1,1)-mu(ii,1)).^2)./((nu(ii,1)-2)*sigma2(1,1));
%         L(1,1) = log(gamma((nu(ii,1)+1)/2))-log(gamma(nu(ii,1)/2))-0.5*log(nu(ii,1)-2) - 0.5*log(pi) ...
%             - 0.5*log(sigma2(1,1))-((nu(ii,1)+1)/2)*log(1+y2);
        
        y2 = ((y(1,1)-mu(ii,1)).^2)./sigma2(1,1);
            
        if (k == 5)
%             pdf(1,1) = duvt_garch(y(1,1), mu(ii,1), rho(ii,1)*sigma2(1,1), nu(ii,1), GamMat);
%             pdf(1,1) = log(pdf(1,1));
            y2 = y2./(nu(ii,1)-2);
            L(1,1) = log(gamma((nu(ii,1)+1)/2)) - log(gamma(nu(ii,1)/2))...
                - 0.5*(log(nu(ii,1)-2) + log(pi) + log(sigma2(1,1)) + (nu(ii,1)+1)*log(1+y2));
        else
            L(1,1) = -0.5*(log(2*pi) + log(sigma2(1,1)) + y2);
        end
            
        for jj = 2:T
            if (k == 5) 
                inv_fisher = nu_con2(ii,1).*(sigma2(jj-1,1).^2)./(fn_chain_rule(sigma2(jj-1,1).^2));             
                w = (nu(ii,1)+1)./(nu(ii,1) - 2 + ((y(jj-1,1)-mu(ii,1)).^2)./sigma2(jj-1,1));
                score =  (w.*((y(jj-1,1)-mu(ii,1)).^2)./sigma2(jj-1,1) - 1)./(2*sigma2(jj-1,1));
            else
                inv_fisher = 2*(sigma2(jj-1,1).^2)./(fn_chain_rule(sigma2(jj-1,1).^2));                 
                score = (((y(jj-1,1)-mu(ii,1)).^2)./sigma2(jj-1,1) - 1)./(2*sigma2(jj-1,1));               
            end
            S = fn_scale(inv_fisher);
            score = fn_chain_rule(sigma2(jj-1,1)).*score;
            scaled_score = S.*score;

            f(jj,1) = omega(ii,1) + A(ii,1)*scaled_score + B(ii,1)*f(jj-1,1);
            sigma2(jj,1) = fn_link(f(jj,1));
            
            y2 = ((y(jj,1)-mu(ii,1)).^2)./sigma2(jj,1);

            if (k == 5)
%                 pdf(jj,1) = duvt_garch(y(jj,1), mu(ii,1), rho(ii,1)*sigma2(jj,1), nu(ii,1), GamMat);
%                 pdf(jj,1) = log(pdf(jj,1));
                y2 = y2./(nu(ii,1)-2);
                L(jj,1) = log(gamma((nu(ii,1)+1)/2)) - log(gamma(nu(ii,1)/2))...
                    - 0.5*(log(nu(ii,1)-2) + log(pi) + log(sigma2(jj,1)) + (nu(ii,1)+1)*log(1+y2));
            else
                L(jj,1) = -0.5*(log(2*pi) + log(sigma2(jj,1)) + y2);
            end
        end
        LL(ii,1) = -sum(L)/T;
%         if (k == 5)
%             d(ii,1) = - sum(pdf)/T;
%             fprintf('d(%i,1) = %8.6f, LL(%i,1) = %8.6f, DIFF = %8.6f\n', ii, d(ii,1), ii, LL(ii,1), d(ii,1)-LL(ii,1));
%         end
    end
end 