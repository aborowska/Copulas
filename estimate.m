function [theta, hessian, hessian_tr, signal_smooth] = estimate(kernel,theta_init,fn_trans_param,fn_jacobian,options) 
    theta_init_trans = fn_trans_param(theta_init, 'opt');
%     NAIS_max = @(par_SV_trans) NAIS_loglik(par_SV_trans, par_NAIS, y, S, cont, RND);
% function [lnL_hat, theta_smooth] = NAIS_loglik(par_SV_trans, par_NAIS, y, S, cont, RND)
%     par_SV =  transform_param_SV(par_SV_trans, [cont.data_on,'_back']);
    
    [theta_trans,~,~,~,~, hessian]= fminunc(kernel, theta_init_trans, options);

    jaco_inv = fn_jacobian(theta_trans);
    hessian_tr = jaco_inv*hessian*jaco_inv;   
    
    if (nargout == 4)
        [~, signal_smooth] = kernel(theta_trans);
    end
    
    theta = fn_trans_param(theta_trans,'back');
end