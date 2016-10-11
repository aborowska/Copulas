function d = posterior_gas(theta, y, L)
    % theta is Nx4, matrix of draws
    [N ,~] = size(theta);
    mu = theta(:,1);
    omega = theta(:,2);
    A = theta(:,3);
    B = theta(:,4); 

    prior = prior_gas(N, omega, B);
    
    T = size(y,1);
     
    d = -Inf*ones(N,1);
       
    for ii = 1:N
        if mod(ii,1000) == 0
            fprintf('posterior ii = %d\n',ii);
        end
        
        f = zeros(T,1);        
        if (prior(ii,1)) % when all the parameter constraints are satisfied
            f(1,1) = omega(ii,1)/(1-B(ii,1)); % unconditional variance to initialize h_1            
            for jj = 2:T
                scaled_score = (y(jj-1,1)-mu(ii,1)).^2 - f(jj-1,1); 
                f(jj,1) = omega(ii,1) + A(ii,1)*scaled_score + B(ii,1)*f(jj-1,1);                       
            end
            pdf = -0.5*(log(2*pi) + log(f) + ((y - mu(ii,1)).^2)./f);
            d(ii,1) = sum(pdf) + prior(ii,2); 
        end
    end
    
    if (~L)
        d = exp(d);
    end
end


function R = prior_gas(N, omega, B)
    % uniform priors 
    
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val an the corresponding point
    
    c1 = (omega > 0);
    c3 = ((B >= 0) & (B < 1));
    
    r1 = (c1 & c3);
    r2 = -Inf*ones(N,1);
    r2(r1==true) = 1; 

    R = [r1, r2];
end
