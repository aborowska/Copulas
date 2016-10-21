addpath(genpath('include/'));
addpath(genpath('KLS/'));

model = 't_gas';

options = optimset('display','iter','TolFun',1e-5,'LargeScale','off','TolX',1e-5,'HessUpdate','bfgs','FinDiffType','central',...
     'maxiter',5000,'MaxFunEvals',5000); %display iter
plot_on = false;

%% EMPIRICAL
    % Given the estimation results for the marginals 
    y = load('ibm_ccola_rets.txt');
    T = size(y,1);
    ibm = y(:,1);
    ccola = y(:,2);

    load('Results/ibm_ccola_results.mat');

    if link  
        fn_link = @(xx) xx;
    else
        fn_link = @(xx) exp(xx);    
    end

    % Standardized residuals, ~t(0,1,nu) (zero mean, unit variance)
    z_ibm = (ibm'-mu_ibm(1,1))./sqrt(fn_link(f_ibm)*(mu_ibm(1,5)-2)/mu_ibm(1,5));
    z_ccola = (ccola'-mu_ccola(1,1))./sqrt(fn_link(f_ccola)*(mu_ccola(1,5)-2)/mu_ccola(1,5));

    if plot_on
        hold on  
        plot(z_ibm,'b') 
        plot(z_ccola,'r')
        hold off
    end
    % cdf transforms of the residuals
    u_ibm = tcdf(z_ibm, mu_ibm(1,5));
    u_ccola = tcdf(z_ccola, mu_ccola(1,5));
    if plot_on
        hold on  
        plot(u_ibm,'b') 
        plot(u_ccola,'r')
        hold off
    end
    
    u = [u_ibm; u_ccola]';
    % mu_init = [0.02, 0.10, 0.98];       
    % [f, rho] = volatility_copula_gas(mu_init, u);

    if plot_on
        hold on  
        % plot(f_ibm,'b') 
        % plot(f_ccola,'r')
        plot(rho,'k') 
        plot(rho*0,'k:')
        hold off
    end


%% SIMULATION


%% Normal Copula estimation

fn_trans_param = @(xx, mm) transform_param_gas_copula(xx, mm);
fn_jacobian = @(xx) jacobian_gas_copula(xx);
link = 1;  % 1: 2*(logsig-0.5) (KLS); 0: inv Fisher transform (Hafner & Manner 2012);
scale = 0; % 1 - inv fisher; 0 - sqrt inv fisher

mu_init = [0.02, 0.10, 0.98];       
kernel_init = @(xx) loglik_copula_gas(fn_trans_param(xx,'back'), u);
tic
[mu_copula, ~, Hessian, signal_copula] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
toc %Elapsed time is 18.795089 seconds.

Sigma_copula = inv(T*Hessian);
[f_copula, rho_copula] = volatility_copula_gas(mu_copula, u);


kernel_init = @(xx) loglik_gen_copula_gas_mex(fn_trans_param(xx,'back'), u, link, scale);
tic
[mu_copula_mex, ~, Hessian_mex, signal_copula_mex] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
toc %Elapsed time is 0.388792 seconds.

if plot_on
    hold on  
    % plot(f_ibm,'b') 
    % plot(f_ccola,'r')
    plot(rho,'b')
    plot(rho_copula,'k') 
    plot(rho*0,'k:')
    hold off
end



%% Student's t Copula estimation
 
mu_init = [0.02, 0.10, 0.98, 8];     
fn_trans_param = @(xx, mm) transform_param_gas_copula(xx, mm);
fn_jacobian = @(xx) jacobian_gas_copula(xx);
kernel_init = @(xx) loglik_copula_t_gas(fn_trans_param(xx,'back'), u);
tic
[mu_copula_t, ~, Hessian, signal_copula_t] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
toc %Elapsed time is 7.022134 seconds.

Sigma_copula_t = inv(T*Hessian);

link = 1; % 1: 2*(logsig-0.5) (KLS); 0: inv Fisher transform (Hafner & Manner 2012);
scale = 0; % 1 - inv fisher; 0 - sqrt inv fisher
kernel_init = @(xx) loglik_gen_copula_gas_mex(fn_trans_param(xx,'back'), u, link, scale);
tic
[mu_copula_t_mex, ~, Hessian_mex, signal_copula_t_mex] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
toc %Elapsed time is 5.247818 seconds.


if plot_on
    hold on  
    % plot(f_ibm,'b') 
    % plot(f_ccola,'r')
    plot(rho_copula,'b') 
    plot(rho_copula_t,'k') 
    plot(rho*0,'k:')
    hold off
end


%% Different link functions:
link = 0; % 1: 2*(logsig-0.5) (KLS); 0: inv Fisher transform (Hafner & Manner 2012);
scale = 0; % 1 - inv fisher; 0 - sqrt inv fisher
kernel_init = @(xx) loglik_gen_copula_gas_mex(fn_trans_param(xx,'back'), u, link, scale);
mu_init = [0.02, 0.10, 0.98];     

tic
[mu_copula_mex2, ~, Hessian_mex2, signal_copula_mex2] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
toc

mu_init = [0.02, 0.10, 0.98, 8];     

tic
[mu_copula_t_mex2, ~, Hessian_mex2, signal_copula_t_mex2] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
toc