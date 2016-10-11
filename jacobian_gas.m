function jaco_inv = jacobian_gas(param_trans, link)
    d = size(param_trans,2);
    
    if (d == 4)
        if link
            jaco_inv = diag([1, 1/exp(param_trans(1,2)), 1, (1+exp(param_trans(1,4)))*(1+exp(-param_trans(1,4)))]);
        else
            jaco_inv = diag([1, 1, 1, (1+exp(param_trans(1,4)))*(1+exp(-param_trans(1,4)))]);
        end
    else
        if link
            jaco_inv = diag([1, 1/exp(param_trans(1,2)), 1, (1+exp(param_trans(1,4)))*(1+exp(-param_trans(1,4))), 1/(exp(param_trans(1,5)))]);
        else
            jaco_inv = diag([1, 1, 1, (1+exp(param_trans(1,4)))*(1+exp(-param_trans(1,4))), 1/(exp(param_trans(1,5)))]);            
        end
    end 
end
    