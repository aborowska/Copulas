clear
clc

S = RandStream('mt19937ar','seed',sum(67*clock));
RandStream.setDefaultStream(S);

%% Settings
filename='results_st_098_n500.mat';
replications=1000; % series realisations
n=2500; % time series length
window=2500; % out of sample period length


%% State Space Model
model= SsfSet('model','SC','p',1,'ARMA',[1 0 1]);
model.psi=0; % number of static parameters
PARAM=[1 0.98 0.01]'; % state space parameters
[c,Z,~,d,T,R,Q,a1,P1]=GetSsf(PARAM,model); % state space
func=(@(x) 2*((1./(1+exp(-x)))-0.5)); % function handle for the parameter
func1=(@(x) 2*((1./(1+exp(-x)))-0.5)); % function handle for the prediction
func2=(@(x) x); 
func3=(@(x) x); % function handle for the prediction

[NUM,Weights]=GHNodes(20); % nodes for the Gauss-Hermite quadrature

%% Starting values for the GAS/GAM models
initial2=[PARAM(1)*(1-0.98) 0.01 0.98]';
initial3=[PARAM(1)*(1-0.98) 0.01 0.98]';
initial4=[0.45*(1-0.98) 0.01 0.98]';
initial5=[0.45*(1-0.98) 0.01 0.98]';
initial6=[0.45*(1-0.98) 0.01 0.90]';

i=1;
counter=0;
errorcounter=0;
observation=zeros(replications*window,1);
target=zeros(replications*window,1);
pred1=zeros(replications*window,6);
pred2=zeros(replications*window,6);
pred3=zeros(replications*window,6);
PARAM_table=zeros(size(PARAM,1),replications,6);
while i<=replications
    disp(i);
    
    %% Main part
    try
    %% State space model simulation
    [y,~,lambda]=simulate(c,Z,d,T,R,Q,a1,P1,n+window);
    y=norminv(y);
    
    %% NAIS estimation
    RND1=normrnd(0,1,100,n);
    PARAM1(:,1)=IsEstim_NAIS(y(:,1:n),model,PARAM,NUM,Weights,RND1(:,1:n));
    PARAM1(:,2)=IsEstim_NAIS(y(:,n-1000:n),model,PARAM1(:,1),NUM,Weights,RND1(:,n-1000:n));
    PARAM1(:,3)=IsEstim_NAIS(y(:,n-500:n),model,PARAM1(:,2),NUM,Weights,RND1(:,n-500:n));
    
    %% GAS/GAM estimations
    PARAM2(:,1)=estimate_gas1(y(:,1:n),initial2);
    PARAM3(:,1)=estimate_acm1(y(:,1:n),initial6);

    PARAM2(:,2)=estimate_gas1(y(:,n-1000:n),PARAM2(:,1));
    PARAM3(:,2)=estimate_acm1(y(:,n-1000:n),PARAM3(:,1));
    
    PARAM2(:,3)=estimate_gas1(y(:,n-500:n),PARAM2(:,2));
    PARAM3(:,3)=estimate_acm1(y(:,n-500:n),PARAM3(:,2));
  
    %% GAS/GAM forecasts
    [~,pred1(counter+1:counter+window,3)]=recursion_gas1(y,PARAM2(:,1),window);
    [~,pred1(counter+1:counter+window,4)]=recursion_acm1(y,PARAM3(:,1),window);
    
    [~,pred2(counter+1:counter+window,3)]=recursion_gas1(y,PARAM2(:,2),window);
    [~,pred2(counter+1:counter+window,4)]=recursion_acm1(y,PARAM3(:,2),window);
  
    [~,pred3(counter+1:counter+window,3)]=recursion_gas1(y,PARAM2(:,3),window);
    [~,pred3(counter+1:counter+window,4)]=recursion_acm1(y,PARAM3(:,3),window);
    
  
     %% Targets
     observation(counter+1:counter+window,1)=y(1,n+1:n+window)';
     target(counter+1:counter+window,1)=lambda(1,n+1:n+window)';
     
     %% IS starting values
     [mu1,V1]=IsStart_NAIS(y,model,PARAM,NUM,Weights);
     [mu2,V2]=IsStart_NAIS(y,model,PARAM1(:,1),NUM,Weights);
     
     %% Forecasting
     for j=1:window
     counter=counter+1;
     RND2=normrnd(0,1,100,250+1);
     
     [pred1(counter,1),~]=IsPred_NAIS_signal(y(:,n+j-250:n+j-1),model,PARAM,NUM,Weights,RND2,func1,func2,func3,mu1(:,n+j-250:n+j-1),V1(:,n+j-250:n+j-1));
     [pred1(counter,2),~]=IsPred_NAIS_signal(y(:,n+j-250:n+j-1),model,PARAM1(:,1),NUM,Weights,RND2,func1,func2,func3,mu2(:,n+j-250:n+j-1),V2(:,n+j-250:n+j-1));
     [pred2(counter,2),~]=IsPred_NAIS_signal(y(:,n+j-250:n+j-1),model,PARAM1(:,2),NUM,Weights,RND2,func1,func2,func3,mu2(:,n+j-250:n+j-1),V2(:,n+j-250:n+j-1));
     [pred3(counter,2),~]=IsPred_NAIS_signal(y(:,n+j-250:n+j-1),model,PARAM1(:,3),NUM,Weights,RND2,func1,func2,func3,mu2(:,n+j-250:n+j-1),V2(:,n+j-250:n+j-1));
     end  
     i=i+1;
     
     if rem(i,50)==0
     save(filename);
     end 
     
     if rem(i,50)==0
     save(filename);
     end 
     
    catch exception
    disp('There was an error');
    disp(exception.stack(1,1).name)
    disp(exception.stack(1,1).line)
    disp(exception.identifier)
    disp(exception.message)
    errorcounter=errorcounter+1;
    if errorcounter==1000
       error('Too many errors, stopping.');
    end
    end
    
end

clear RND1 RND2 NUM Weights i j y y1 y2 ytilde1 ytilde2 H1 H2 lambda
clear initial2 initial3 initial4 initial5 initial6 PARAM1 PARAM2 PARAM2 PARAM3 PARAM4 PARAM5 PARAM6
save(filename);

pred2=pred2(:,2:end);
pred3=pred3(:,2:end);

e=repmat(target(:,1),1,size(pred1,2))-pred1;
disp([mean(e);var(e);mean(e.^2);]')

e=repmat(target(:,1),1,size(pred2,2))-pred2;
disp([mean(e);var(e);mean(e.^2);]')

e=repmat(target(:,1),1,size(pred3,2))-pred3;
disp([mean(e);var(e);mean(e.^2);]')



