function [y_hp, eps_hp, h] = predict_t_garch_noS(theta, y_T, h_T, hp, eps)
    [N ,~] = size(theta);
    omega = theta(:,1);
    alpha = theta(:,2);
    beta = theta(:,3);
    mu = theta(:,4);
    nu = theta(:,5);
      
    rho = (nu-2)./nu;
    
    if (nargin == 4)
%         fprintf('hp = %i \n',hp);
        eps_hp = trnd(repmat(nu,1,hp));
    else %(with given eps)
        eps_hp = eps;
    end
    
    y_hp = zeros(N,hp+1);    
    y_hp(:,1) = y_T.*ones(N,1);

    h = zeros(N,hp+1); 
    h(:,1) = h_T;
    
    for jj = 2:(hp+1)
         h(:,jj) = omega(:,1) + alpha(:,1).*(y_hp(:,jj-1)-mu(:,1)).^2 + beta(:,1).*h(:,jj-1);  
         y_hp(:,jj) = mu(:,1) + sqrt(rho(:,1).*h(:,jj)).*eps_hp(:,jj-1);
    end
    y_hp = y_hp(:,2:hp+1);
    h = h(:,2:hp+1);
%     y_hp = predict_t_garch_noS_mex(theta, y_T, h_T, eps_hp);
end