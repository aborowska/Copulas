function [u,theta,rho]=simulate(c,Z,d,T,R,Q,a1,P1,n)
    m=size(T,1);
    r=size(Q,1);
    theta=SsfStateSim(c,Z,d,T,R,Q,a1,P1,n);

    rho=2*(logsig(theta)-0.5);

    u1=zeros(1,n);
    u2=zeros(1,n);

    for i=1:n
        U=copularnd('Gaussian',rho(1,i),1);
        u1(1,i)=U(1,1);
        u2(1,i)=U(1,2);
    end

    u=[u1;u2];
end