%% Gaussian
% With anonymous functions
mu_init = [0.02, 0.10, 0.98];       
kernel_init = @(xx) loglik_copula_gas(fn_trans_param(xx,'back'), u);
tic
% profile on
[mu_copula, ~, Hessian, signal_copula] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
% profile off
% profile viewer
toc %Elapsed time is 18.795089 seconds.



% With subfunction
mu_init = [0.02, 0.10, 0.98];       
kernel_init = @(xx) loglik_copula_gas2(fn_trans_param(xx,'back'), u);
% tic
profile on
[mu_copula2, ~, Hessian2, signal_copula2] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
profile off
profile viewer
% toc %Elapsed time is 8.523910 seconds.

%% Student's t 
% With anonymous functions
mu_init = [0.02, 0.10, 0.98, 8];       
kernel_init = @(xx) loglik_copula_t_gas(fn_trans_param(xx,'back'), u);
% tic
profile on
[mu_copula, ~, Hessian, signal_copula] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
profile off
profile viewer
% toc %Elapsed time is 7.028290 sseconds.



% With subfunction
mu_init = [0.02, 0.10, 0.98, 8];       
kernel_init = @(xx) loglik_copula_t_gas2(fn_trans_param(xx,'back'), u);
% tic
% profile on
[mu_copula2, ~, Hessian2, signal_copula2] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
% profile off
% profile viewer
toc %Elapsed time is 23.376660 seconds.

