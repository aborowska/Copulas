function h_T = volatility_t_garch_noS(theta, data, S)
% function h_T = volatility_t_garch(theta, data, S)

    [N ,~] = size(theta);
    omega = theta(:,1);
    alpha = theta(:,2);
    beta = theta(:,3);
    mu = theta(:,4);

    T = size(data,1);
%     ind = 2:T;
    data = data';
    
    h = zeros(N,T); 
    h(:,1) = S*ones(N,1);
    
%     h(:,ind) = repmat(omega(:,1),1,T-1) + repmat(alpha(:,1),1,T-1).*(repmat(data(1,ind-1),N,1)-repmat(mu(:,1),1,T-1)).^2;
    for jj = 2:T
        h(:,jj) = omega(:,1) + alpha(:,1).*(data(1,jj-1) - mu(:,1)).^2 + beta(:,1).*h(:,jj-1) ;
    end
    h_T = h(:,T);
end