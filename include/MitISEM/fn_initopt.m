function [mu, Sigma, val] = fn_initopt(kernel, mu0, fn_delta) % , options)
    d = length(mu0);
%     options = optimset('TolX', 0.0001, 'Display', 'iter', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');

    [mu,val,~,~,~,hessian] = fminunc(kernel, mu0); % ,options);
%     options = optimset('Display','iter');
%     x = fminsearch(kernel_init,mu_init, options);
%     x = fminsearch(kernel_init,mu_hl), options);
%     [x,~,~,~,~,hessian] = fminunc(kernel,mu0,options)
    if nargin > 2
        hessian = fn_delta(mu, hessian);
    end
    try
        [~, T] = kernel(mu0);
        Sigma = inv(T*hessian);
    catch
        Sigma = inv(hessian);
    end
    Sigma = reshape(Sigma,1,d^2);  
end