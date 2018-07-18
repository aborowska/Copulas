function [u_sim, f_true, rho_true] = simulate_copula_gas(mu_true, T, N, link)
    d = size(mu_true,2); 

    omega = mu_true(:,1);
    A = mu_true(:,2);
    B = mu_true(:,3);
    if (d == 4)
        nu = mu_true(:,4);
    end
    if link
        transf = @(aa) (1 - exp(-aa))./(1 + exp(-aa));
    else
        transf = @(aa) (exp(2*aa)-1)./(exp(2*aa)+1);        
    end 
    f_true = zeros(N,T);
    rho_true = zeros(N,T);    
    u_sim = zeros(T,2,N);   
    

    f(:,1) = omega./(1-B);
    rho_true(:,1) = transf(f(:,1));
    if (d == 3)
        scoref = @(z1, z2, r) ((1 + r.^2).*(z1.*z2 - r) - r.*(z1.^2 + z2.^2 - 2))./((1 - r.^2).^2);
        inff = @(r) (1 + r.^2)./((1 - r.^2).^2);

        for jj = 2:T
            u_sim(jj-1,:,ii) = bi_copula_n_rnd(rho_true(jj-1,:));
            z = norminv(u_sim(jj-1,:,ii));
            s = scoref(z(1,1), z(1,2), rho_true(:,jj-1));
            scscore = s./sqrt(inff(rho_true(:,jj-1)));
            f_true(:,jj) = omega + A.*scscore + B.*f_true(:,jj-1);
            rho_true(:,jj) = transf(f_true(:,jj));        
        end
    else
        scoref = @(z1, z2, r, w) ((1 + r.^2).*(w.*z1.*z2 - r) - r.*(w.*z1.^2 + w.*z2.^2 - 2))./((1 - r.^2).^2);
        inff = @(r, v) (v + 2 + v.* r.^2)./((v+4).*(1 - r.^2).^2);
        wf = @(z1, z2, r, v) (v + 2)./(v + (z1.^2 + z2.^2 - 2*r.*z1.*z2)./(1 - r.^2)) ;

        for ii = 1:N
            z = tinv(y,nu(ii,1));

            for jj = 2:T
                            u_sim(:,:,jj) = bi_copula_n_rnd(rho_true(jj,:));

                u_sim(jj-1,:) = bi_copula_t_rnd(rho_true(jj-1,:), nu);
                z = tinv(u_sim(jj-1,:),nu);

                w = wf(z(1,1), z(1,2), rho_true(ii,jj-1), nu(ii,1));
                s = scoref(z(1,1), z(1,2), rho_true(ii,jj-1), w);
                scscore = s./sqrt(inff(rho_true(ii,jj-1), nu(ii,1)));

                f_true(ii,jj) = omega(ii,1) + A(ii,1).*scscore + B(ii,1).*f_true(ii,jj-1);
                rho_true(ii,jj) = transf(f_true(ii,jj));        
            end
        end     
    end        
 
    if (N == 1)
        u_sim = squeeze(u_sim);
    end 
end
