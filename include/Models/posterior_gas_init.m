function d = posterior_gas_init(theta, y)
    mu = theta(:,1);
    omega = exp(theta(:,2));
    A = theta(:,3);
    B = logsig(theta(:,4)); 
    
    T = size(y,1);
             
    f = zeros(T,1);
    f(1,1) = omega/(1-B); % unconditional variance to initialize h_1            
    for t = 2:T
        scaled_score = (y(t-1,1)-mu).^2 - f(t-1,1); 
        f(t,1) = omega + A*scaled_score + B*f(t-1,1);                       
    end
    pdf = -0.5*(log(2*pi) + log(f) + ((y - mu).^2)./f);
    d = sum(pdf)/T; 
end