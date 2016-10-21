function jaco_inv = jacobian_garch(param_trans)
    d = size(param_trans,2);
    
    if (d == 4)
        jaco_inv = diag([1/exp(param_trans(1,1)), (1+exp(param_trans(1,2)))*(1+exp(-param_trans(1,2))),...
            (1+exp(param_trans(1,3)))*(1+exp(-param_trans(1,3))), 1]);       
    else
        jaco_inv = diag([1/exp(param_trans(1,1)), (1+exp(param_trans(1,2)))*(1+exp(-param_trans(1,2))), ...
            (1+exp(param_trans(1,3)))*(1+exp(-param_trans(1,3))), 1, 1/(exp(param_trans(1,5)))]);
    end 
end
    