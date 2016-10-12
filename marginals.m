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

%% ibm
link =  0; % 1 - linear link; 0 - exp link
scale = 0; % 1 - inv fisher; 0 - sqrt inv fisher
% Note: (1,1) was discussed in the original GAS paper (Creal et al. 2013)
% Note: (0,0) cancels out so the updating step has a constant unit variance
% and is invariant under any nondegenerate parameter transformation 

% theta = [mu, omega, A, B, nu] 
theta_init = [0.07, 0.01, 0.1, 0.98, 15]; 
d = size(theta_init,2);
  
% kernel_init = @(xx) posterior_t_gas_hyper_init_mex(xx, ibm, hyper, GamMat);
% [theta_init_ibm, Sigma_ibm] = fn_initopt(kernel_init, theta_init);
[lnk, f_ibm] = posterior_gen_gas_mex(theta_init, ibm, hyper, link, scale, GamMat);


% Real MLE
fn_trans_param = @(xx,mm) transform_param_gas(xx, mm, link);
kernel_init = @(xx) loglik_gen_gas(fn_trans_param(xx,'back'), ibm, link, scale, GamMat); 
% kernel_init = @(xx) loglik_t_gas_hyper_init_mex(transform_param_gas(xx,'back'), ibm, hyper, GamMat);
options = optimset('display','off','TolFun',1e-5,'LargeScale','off','TolX',1e-5,'HessUpdate','bfgs','FinDiffType','central',...
     'maxiter',5000,'MaxFunEvals',5000); %display iter
fn_jacobian = @(xx) jacobian_gas(xx, link);
[theta_ibm, ~, Hessian, signal_smooth] = estimate(kernel_init, theta_init, fn_trans_param, fn_jacobian, options);
Sigma_ibm = inv(T*Hessian);

mit_ibm = struct('mu',theta_ibm,'Sigma',reshape(Sigma_ibm,1,d^2),'p',1,'df',5);
% kernel = @(xx) posterior_t_gas_hyper_mex(xx, ibm, hyper, GamMat);
kernel = @(xx) posterior_gen_gas_mex(xx, ibm, hyper, link, scale, GamMat);
[theta_ibm, accept_ibm] = Mit_MH(M+BurnIn, kernel, mit_ibm, GamMat);
theta_ibm = theta_ibm(BurnIn+1:M+BurnIn,:);
% accept_ibm = 0.503181818181818


%% DEBUGGING
% PARAM = mu_init(2:5);
% ibm_kls = ibm'-mean(ibm);
PARAM = theta_init;
ibm_kls = ibm';
[PARAM,f,theta]=estimate_vol_t_gas(ibm_kls,PARAM);

[LL1, f1] = loglik_t_gas(PARAM, ibm, link, scale, GamMat);
[LL2, f2] = loglik_gen_gas(PARAM, ibm, link, scale, GamMat);
[LL3, f3] = posterior_gen_gas_mex(PARAM, ibm, hyper, link, scale, GamMat);



%% ccola
kernel_init = @(xx) posterior_t_gas_hyper_init_mex(xx, ccola, hyper, GamMat);
mu_init = [0.07, 0.02, 0.1, 0.98, 5];
mu_init= [ 0.0524    0.0892    0.0550    0.9751    4.8330];

[mu_ccola, ~] = fn_initopt(kernel_init, mu_init);
[mu_ccola, Sigma_ccola] = fn_initopt(kernel_init, mu_ccola);

mit_ccola = struct('mu',mu_ccola,'Sigma',Sigma_ccola,'p',1,'df',5);
kernel = @(xx) posterior_t_gas_hyper_mex(xx, ccola, hyper, GamMat);
[theta_ccola, accept_ccola] = Mit_MH(M+BurnIn, kernel, mit_ccola, GamMat);
theta_ccola = theta_ccola(BurnIn+1:M+BurnIn,:);
  