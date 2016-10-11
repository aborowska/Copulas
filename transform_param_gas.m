function param_trans = transform_param_gas(param, mode, link)
    param_trans = param;
    d = size(param,2);
    
    switch mode
%         case 'sim_opt' % transformation for optimization --> unbounded
%             param(1,2) = param(1,2)/2 + 0.5;
%             param_trans(1,2) = log(param(1,2)/(1-param(1,2)));
%             param_trans(1,3) = log(param(1,3));
%             if (d == 4)
%                 param_trans(1,4) = log(param(1,4) - 2);
%             end
%             
%         case 'sim_back' % if mode == 'back' (transform back)
%             param_trans(1,2) = exp(param(1,2))/(1+exp(param(1,2)));
%             param_trans(1,2) = 2*(param_trans(1,2) - 0.5);
%             param_trans(1,3) = exp(param(1,3));
%             if (d == 4)
%                 param_trans(1,4) = 2 + exp(param(1,4));
%             end
            
%         case 'est_opt' % transformation for optimization --> unbounded
        case 'opt' % transformation for optimization --> unbounded
            if link 
                param_trans(1,2) = log(param(1,2));
            end
            param_trans(1,4) = log(param(1,4)/(1-param(1,4)));
            if (d == 5)
                param_trans(1,5) = log(param(1,5) - 2);
            end
            
%         case 'est_back'% if mode == 'back' (transform back)
        case 'back'% if mode == 'back' (transform back)
            if link 
                param_trans(1,2) = exp(param(1,2));
            end
            param_trans(1,4) = exp(param(1,4))/(1+exp(param(1,4)));
            param_trans(1,4) = min(0.998, param_trans(1,4));
            if (d == 5)
                param_trans(1,5) = 2 + exp(param(1,5));
                param_trans(1,5) = min(200, param_trans(1,5));
            end
    end
end