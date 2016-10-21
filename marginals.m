clear all
close all
addpath(genpath('include/'));
addpath(genpath('KLS/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

M = 10000;
BurnIn = 1000;

y = load('ibm_ccola_rets.txt');
T = size(y,1);
ibm = y(:,1);
ccola = y(:,2);

model = 't_gas';

hyper = 0.01; % for posterior: prior hyperparameter for nu, i.e. df 

options = optimset('display','off','TolFun',1e-5,'LargeScale','off',...
    'TolX',1e-5,'HessUpdate','bfgs','FinDiffType','central',...
     'maxiter',5000,'MaxFunEvals',5000); %display iter
 
%% ibm
link =  0; % 1 - linear link; 0 - exp link
scale = 0; % 1 - inv fisher; 0 - sqrt inv fisher
% Note: (1,1) was discussed in the original GAS paper (Creal et al. 2013)
% Note: (0,0) cancels out so the updating step has a constant unit variance
% and is invariant under any nondegenerate parameter transformation 
fn_trans_param = @(xx,mm) transform_param_gas(xx, mm, link);
fn_jacobian = @(xx) jacobian_gas(xx, link);
if link  
    fn_link = @(xx) xx;
else
    fn_link = @(xx) exp(xx);    
end

% theta = [mu, omega, A, B, nu] 
mu_init = [0.07, 0.01, 0.1, 0.98, 15]; 
d = size(mu_init,2);
  
% kernel_init = @(xx) posterior_t_gas_hyper_init_mex(xx, ibm, hyper, GamMat);
% [theta_init_ibm, Sigma_ibm] = fn_initopt(kernel_init, theta_init);
% [lnk, f_ibm] = posterior_gen_gas_mex(theta_init, ibm, hyper, link, scale, GamMat);


% MLE
kernel_init = @(xx) loglik_gen_gas_mex(fn_trans_param(xx,'back'), ibm, link, scale); 
% kernel_init = @(xx) loglik_gen_gas(fn_trans_param(xx,'back'), ibm, link, scale, GamMat); 
% kernel_init = @(xx) loglik_t_gas_hyper_init_mex(transform_param_gas(xx,'back'), ibm, hyper, GamMat);
tic
[mu_ibm, ~, Hessian, signal_ibm] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
toc
Sigma_ibm = inv(T*Hessian);

mit_ibm = struct('mu',mu_ibm,'Sigma',reshape(Sigma_ibm,1,d^2),'p',1,'df',5);
% kernel = @(xx) posterior_t_gas_hyper_mex(xx, ibm, hyper, GamMat);
kernel = @(xx) posterior_gen_gas_mex(xx, ibm, hyper, link, scale, GamMat);
[theta_ibm, accept_ibm] = Mit_MH(M+BurnIn, kernel, mit_ibm, GamMat);
theta_ibm = theta_ibm(BurnIn+1:M+BurnIn,:);
% accept_ibm = 0.503181818181818


%% DEBUGGING
% PARAM = mu_init(2:5);
% ibm_kls = ibm'-mean(ibm);
PARAM = mu_init;
ibm_kls = ibm';
[PARAM,f,sigma2]=estimate_vol_t_gas(ibm_kls,PARAM);

PARAM = [mu_ibm;mu_init];
[LL0,f0]=loglik_vol_t_gas(ibm_kls,PARAM);

[LL1, f1] = loglik_t_gas(PARAM, ibm, link, scale, GamMat);
[LL2, f2] = loglik_gen_gas(PARAM, ibm, link, scale, GamMat);
[LL3, f3] = loglik_gen_gas_mex(PARAM, ibm, link, scale);

[LL4, f4] = posterior_gen_gas_mex(PARAM, ibm, hyper, link, scale, GamMat);

hold on
plot(f0-f3,'k')
plot(f3,'r')
hold off


% Speed comparison:
fn_trans_param = @(xx,mm) transform_param_gas(xx, mm, link);
kernel_init = @(xx) loglik_gen_gas(fn_trans_param(xx,'back'), ibm, link, scale, GamMat); 
% kernel_init = @(xx) loglik_t_gas_hyper_init_mex(transform_param_gas(xx,'back'), ibm, hyper, GamMat);
options = optimset('display','off','TolFun',1e-5,'LargeScale','off','TolX',1e-5,'HessUpdate','bfgs','FinDiffType','central',...
     'maxiter',5000,'MaxFunEvals',5000); %display iter
fn_jacobian = @(xx) jacobian_gas(xx, link);
tic
[mu_ibm, ~, Hessian, signal_ibm] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
toc
% Elapsed time is 37.392040 seconds.

% KLS
ibm_kls = ibm';
PARAM = mu_init;
tic
[PARAM,f]=estimate_vol_t_gas(ibm_kls,PARAM,options);
toc
% Elapsed time is 6.338096 seconds.

% My MEXED
kernel_init = @(xx) loglik_gen_gas_mex(fn_trans_param(xx,'back'), ibm, link, scale); 
tic 
[mu_ibm2, ~, Hessian2, signal_ibm2] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
toc
% Elapsed time is 0.509850 seconds.

% RL:
% Elapsed time is 56.003601 seconds.
vp_mle = [ 0.0183, 0.0791, 0.9772, 0.0794,8.5297];
GasVolaUnivMain

%% ccola
kernel_init = @(xx) loglik_gen_gas_mex(fn_trans_param(xx,'back'), ccola, link, scale); 
mu_init = [0.07, 0.01, 0.1, 0.98, 15]; 
% mu_init = [0.07, 0.02, 0.1, 0.98, 5];
% mu_init= [ 0.0524    0.0892    0.0550    0.9751    4.8330];
tic
[mu_ccola, ~, Hessian, signal_ccola] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
toc
Sigma_ccola = inv(T*Hessian); 

mit_ccola = struct('mu',mu_ccola,'Sigma',reshape(Sigma_ccola,1,d^2),'p',1,'df',5);
kernel = @(xx) posterior_gen_gas_mex(xx, ccola, hyper, link, scale, GamMat);
[theta_ccola, accept_ccola] = Mit_MH(M+BurnIn, kernel, mit_ccola, GamMat);
theta_ccola = theta_ccola(BurnIn+1:M+BurnIn,:);
%  accept_ccola = 0.

%% statistics 
% mean(theta_ccola)
% mean(theta_ibm)

[~, f_ibm] = loglik_gen_gas_mex(mu_ibm, ibm, link, scale);
[~, f_ccola] = loglik_gen_gas_mex(mu_ccola, ccola, link, scale);

hold on 
% plot(ibm',':c')
% plot(ccola',':m')
plot(f_ibm,'b')
plot(f_ccola,'r')
hold off
  

save('Results/ibm_ccola_results.mat', 'link', 'scale', 'mu_ibm', 'mu_ccola',...
    'Sigma_ibm', 'Sigma_ccola', 'f_ibm', 'f_ccola');











%% Crisis SnP500 data (PMitISEM project)
snp = csvread('GSPC_ret_updated.csv'); 
snp = snp*100;
            mu_init = [0.07, 0.02, 0.1, 0.98, 6.9];
            %     [val, Td] = kernel_init(mu_init);
            tic
            [mu, Sigma] = fn_initopt(kernel_init, mu_init);           
            mit_direct = struct('mu',mu,'Sigma',Sigma,'p',1,'df',cont_direct.mit.dfnc);
            time_direct(1,1) = toc;   
link = 1;
scale = 1;
mu_init = [0.07, 0.01, 0.1, 0.98, 15]; 
d = size(mu_init,2);
 
[lnk_snp, f_snp] = posterior_gen_gas_mex(mu_init, snp, hyper, link, scale, GamMat);

% Real MLE
fn_trans_param = @(xx,mm) transform_param_gas(xx, mm, link);
kernel_init = @(xx) loglik_gen_gas_mex(fn_trans_param(xx,'back'), snp, link, scale); 
 % kernel_init = @(xx) loglik_t_gas_hyper_init_mex(transform_param_gas(xx,'back'), ibm, hyper, GamMat);
options = optimset('display','off','TolFun',1e-5,'LargeScale','off','TolX',1e-5,'HessUpdate','bfgs','FinDiffType','central',...
     'maxiter',5000,'MaxFunEvals',5000); %display iter
fn_jacobian = @(xx) jacobian_gas(xx, link);
[mu_snp, ~, Hessian, signal_snp] = estimate(kernel_init, mu_init, fn_trans_param, fn_jacobian, options);
Sigma_snp = inv(T*Hessian);