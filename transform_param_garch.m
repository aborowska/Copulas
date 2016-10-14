function param_trans = transform_param_garch(param, mode)
    param_trans = param;
    d = size(param,2);
    
    switch mode            
%         case 'est_opt' % transformation for optimization --> unbounded
        case 'opt' % transformation for optimization --> unbounded
            param_trans(1,1) = log(param(1,1));
            param_trans(1,2) = log(param(1,2)/(1-param(1,2)));
            param_trans(1,3) = log(param(1,3)/(1-param(1,3)));
            if (d == 5)
                param_trans(1,5) = log(param(1,5) - 2);
            end
            
%         case 'est_back'% if mode == 'back' (transform back)
        case 'back'% if mode == 'back' (transform back)
            param_trans(1,1) = exp(param(1,1));
            param_trans(1,2) = exp(param(1,2))/(1+exp(param(1,2)));
            param_trans(1,3) = exp(param(1,3))/(1+exp(param(1,3)));
            param_trans(1,3) = min(0.998, param_trans(1,3));            
            if (d == 5)
                param_trans(1,5) = 2 + exp(param(1,5));
                param_trans(1,5) = min(200, param_trans(1,5));
            end
    end
end