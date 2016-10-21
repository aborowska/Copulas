addpath(genpath('include/'));
addpath(genpath('KLS/'));


y = load('ibm_ccola_rets.txt');
T = size(y,1);
ibm = y(:,1);
ccola = y(:,2);

model = 't_gas';

plot_on = false;

%% Given the estimation results for the marginals  
load('Results/ibm_ccola_results.mat');

% fn_trans_param = @(xx,mm) transform_param_gas(xx, mm, link);
% fn_jacobian = @(xx) jacobian_gas(xx, link);
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
theta = [0.02, 0.10, 0.98];       
[f, rho] = volatility_copula_gas(theta, u);

if plot_on
    hold on  
    % plot(f_ibm,'b') 
    % plot(f_ccola,'r')
    plot(rho,'k') 
    plot(rho*0,'k:')
    hold off
end

%% Normal Copula estimation

fn_trans_param = @(xx, mm) transform_param_gas_copula(xx, mm);
fn_jacobian = @(xx) jacobian_gas_copula(xx);


mu_init = [0.02, 0.10, 0.98];       
kernel_init = @(xx) loglik_copula_gas(fn_trans_param(xx,'back'), u);
options = optimset('display','iter','TolFun',1e-5,'LargeScale','off','TolX',1e-5,'HessUpdate','bfgs','FinDiffType','central',...
     'maxiter',5000,'MaxFunEvals',5000); %display iter
tic
[mu_copula, ~, Hessian, signal_copula] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
toc
Sigma_copula = inv(T*Hessian);

[f_copula, rho_copula] = volatility_copula_gas(mu_copula, u);

kernel_init = @(xx) loglik_gen_copula_gas_mex(fn_trans_param(xx,'back'), u, link, scale);
tic
[mu_copula_mex, ~, Hessian_mex, signal_copula_mex] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
toc

if plot_on
    hold on  
    % plot(f_ibm,'b') 
    % plot(f_ccola,'r')
    plot(rho,'b')
    plot(rho_copula,'k') 
    plot(rho*0,'k:')
    hold off
end

link = 1;  % 1: KLS link; 0: NAIS link;
scale = 0; % 1 - inv fisher; 0 - sqrt inv fisher
mu_init = [0.02, 0.10, 0.98];       
kernel_init = @(xx) loglik_gen_copula_gas_mex(fn_trans_param(xx,'back'), u, link, scale);

[LL1, f1, rho1] = loglik_copula_gas(mu_init, u);
[LL2, f2, rho2] = loglik_gen_copula_gas_mex(mu_init, u, link, scale );



%% Student's t Copula estimation
 
mu_init = [0.02, 0.10, 0.98, 8];     
fn_trans_param = @(xx, mm) transform_param_gas_copula(xx, mm);
fn_jacobian = @(xx) jacobian_gas_copula(xx);
% [LL_copula_t, f_copula_t, rho_copula_t] = loglik_copula_t_gas(mu_init, u);
kernel_init = @(xx) loglik_copula_t_gas(fn_trans_param(xx,'back'), u);
options = optimset('display','iter','TolFun',1e-5,'LargeScale','off','TolX',1e-5,'HessUpdate','bfgs','FinDiffType','central',...
     'maxiter',5000,'MaxFunEvals',5000); %display iter
tic
[mu_copula_t, ~, Hessian, signal_copula_t] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
toc
Sigma_copula_t = inv(T*Hessian);

[~, f_copula_t, rho_copula_t] = loglik_copula_t_gas(mu_copula_t, u);
if plot_on
    hold on  
    % plot(f_ibm,'b') 
    % plot(f_ccola,'r')
    plot(rho_copula,'b') 
    plot(rho_copula_t,'k') 
    plot(rho*0,'k:')
    hold off
end



[LL1, f1, rho1] = loglik_copula_t_gas(mu_init, u);
[LL2, f2, rho2] = loglik_gen_copula_gas_mex(mu_init, u, link, scale );

