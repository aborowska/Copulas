function [LL, f, sigma2] = loglik_t_gas(theta, y, link, scale, GamMat)
    % theta is Nx5, matrix of draws
    [N ,~] = size(theta);
    mu = theta(:,1);
    omega = theta(:,2);
    A = theta(:,3);
    B = theta(:,4); 
    nu = theta(:,5);
    
    rho = (nu-2)./nu;
    
%     nu_con1 = (nu+1)./(nu-2);
    nu_con2 = 2*((nu+3)./nu);

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
     
    d = -Inf*ones(N,1);
    LL = -Inf*ones(N,1);

    for ii = 1:N
        if mod(ii,1000) == 0
            fprintf('loglik ii = %d\n',ii);
        end
        
        f = zeros(T,1);
        sigma2 = zeros(T,1);
        
        pdf = zeros(T,1);
        L = zeros(T,1);

        f(1,1) = omega(ii,1)/(1-B(ii,1)); % unconditional variance to initialize h_1
        sigma2(1,1) = fn_link(f(1,1));

% %         pdf(1,1) = dmvt(y(1,1), mu(ii,1), rho(ii,1)*sigma2(1,1), nu(ii,1), GamMat);
        pdf(1,1) = duvt_garch(y(1,1), mu(ii,1), rho(ii,1)*sigma2(1,1), nu(ii,1), GamMat);
        pdf(1,1) = log(pdf(1,1));
        y2 = ((y(1,1)-mu(ii,1)).^2)./((nu(ii,1)-2)*sigma2(1,1));
        L(1,1) = log(gamma((nu(ii,1)+1)/2))-log(gamma(nu(ii,1)/2))-0.5*log(nu(ii,1)-2) - 0.5*log(pi) ...
            - 0.5*log(sigma2(1,1))-((nu(ii,1)+1)/2)*log(1+y2);
            
        for jj = 2:T
            inv_fisher = nu_con2(ii,1).*(sigma2(jj-1,1).^2)./(fn_chain_rule(sigma2(jj-1,1).^2));             
            scale = fn_scale(inv_fisher);
            w = (nu(ii,1)+1)./(nu(ii,1) - 2 + ((y(jj-1,1)-mu(ii,1)).^2)./sigma2(jj-1,1));
            score =  (w.*((y(jj-1,1)-mu(ii,1)).^2)./sigma2(jj-1,1) - 1)./(2*sigma2(jj-1,1));
            score = fn_chain_rule(sigma2(jj-1,1)).*score;
            scaled_score = scale.*score;

            f(jj,1) = omega(ii,1) + A(ii,1)*scaled_score + B(ii,1)*f(jj-1,1);
            sigma2(jj,1) = fn_link(f(jj,1));

% %                 pdf(jj,1) = dmvt(y(jj,1), mu(ii,1), rho(ii,1)*f(jj,1), nu(ii,1), GamMat);
            pdf(jj,1) = duvt_garch(y(jj,1), mu(ii,1), rho(ii,1)*sigma2(jj,1), nu(ii,1), GamMat);
            pdf(jj,1) = log(pdf(jj,1));
            y2 = ((y(jj,1)-mu(ii,1)).^2)./((nu(ii,1)-2)*sigma2(jj,1));
            L(jj,1) = log(gamma((nu(ii,1)+1)/2))-log(gamma(nu(ii,1)/2))-0.5*log(nu(ii,1)-2)- 0.5*log(pi) ...
                - 0.5*log(sigma2(jj,1))-((nu(ii,1)+1)/2)*log(1+y2);
        end
        d(ii,1) = -sum(pdf)/T;
        LL(ii,1) = -sum(L)/T;
        fprintf('d(%i,1) = %8.6f, LL(%i,1) = %8.6f, DIFF = %8.6f\n', ii, d(ii,1), ii, LL(ii,1), d(ii,1)-LL(ii,1));
    end
end 