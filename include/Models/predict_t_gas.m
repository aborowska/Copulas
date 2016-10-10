function [y_hp, eps_hp, f] = predict_t_gas(theta, y_T, f_T, hp, eps)
    [N ,~] = size(theta);
    mu = theta(:,1);
    omega = theta(:,2);
    A = theta(:,3);
    B = theta(:,4);
    nu = theta(:,5);
      
    rho = (nu-2)./nu;    
    nu_con = (nu+1)./(nu-2);
    A = A.*((nu+3)./nu);
    
    if (nargin == 4)
         eps_hp = trnd(repmat(nu,1,hp));
    else %(with given eps)
        eps_hp = eps;
    end
    
    y_hp = zeros(N,hp+1);    
    y_hp(:,1) = y_T.*ones(N,1);

    f = zeros(N,hp+1); 
    f(:,1) = f_T;
    
    for jj = 2:(hp+1)
         C = 1 + ((y_hp(:,jj-1)-mu).^2)./((nu-2).*f(:,jj-1));              
         f(:,jj) = omega + A.*(nu_con.*((y_hp(:,jj-1)-mu).^2)./C - f(:,jj-1)) + B.*f(:,jj-1);                        
%          f(:,jj) = omega(:,1) + alpha(:,1).*(y_hp(:,jj-1)-mu(:,1)).^2 + beta(:,1).*h(:,jj-1);       
         y_hp(:,jj) = mu(:,1) + sqrt(rho(:,1).*f(:,jj)).*eps_hp(:,jj-1);
    end
    y_hp = y_hp(:,2:hp+1);
    f = f(:,2:hp+1);
end