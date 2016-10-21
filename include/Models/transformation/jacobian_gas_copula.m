function jaco_inv = jacobian_gas_copula(param_trans)
    d = size(param_trans,2);
    
    if (d == 3)
        jaco_inv = diag([1/exp(param_trans(1,1)), 1, ...
            (1+exp(param_trans(1,3)))*(1+exp(-param_trans(1,3)))]);       
    else
        jaco_inv = diag([1/exp(param_trans(1,1)), 1, ...
            (1+exp(param_trans(1,3)))*(1+exp(-param_trans(1,3))), 1/(exp(param_trans(1,4)))]);
    end 
end
    