addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

M = 10000;
BurnIn = 1000;

y = load('ibm_ccola_rets.txt');

ibm = y(:,1);
ccola = y(:,2);

model = 't_gas';
hyper = 0.01;

% ibm
kernel_init = @(xx) posterior_t_gas_hyper_init_mex(xx, ibm, hyper, GamMat);
mu_init = [0.07, 0.02, 0.1, 0.98, 15];
[mu_ibm, Sigma_ibm] = fn_initopt(kernel_init, mu_init);
mit_ibm = struct('mu',mu_ibm,'Sigma',Sigma_ibm,'p',1,'df',10);
kernel = @(xx) posterior_t_gas_hyper_mex(xx, ibm, hyper, GamMat);
[theta_ibm, accept_ibm] = Mit_MH(M+BurnIn, kernel, mit_ibm, GamMat);
theta_ibm = theta_ibm(BurnIn+1:M+BurnIn,:);

% ccola
kernel_init = @(xx) posterior_t_gas_hyper_init_mex(xx, ccola, hyper, GamMat);
mu_init = [0.07, 0.02, 0.1, 0.98, 5];
mu_init= [ 0.0524    0.0892    0.0550    0.9751    4.8330];

[mu_ccola, ~] = fn_initopt(kernel_init, mu_init);
[mu_ccola, Sigma_ccola] = fn_initopt(kernel_init, mu_ccola);

mit_ccola = struct('mu',mu_ccola,'Sigma',Sigma_ccola,'p',1,'df',5);
kernel = @(xx) posterior_t_gas_hyper_mex(xx, ccola, hyper, GamMat);
[theta_ccola, accept_ccola] = Mit_MH(M+BurnIn, kernel, mit_ccola, GamMat);
theta_ccola = theta_ccola(BurnIn+1:M+BurnIn,:);
  