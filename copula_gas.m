 y = load('ibm_ccola_rets.txt');
 theta = [0.02, 0.10, 0.98;
          0.017, 0.08, 0.97];
u = ...; % % transformed series
[f, rho] = volatility_copula_gas(theta, u);

 theta_start = theta(1,:);     
 