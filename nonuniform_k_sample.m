function [S]=nonuniform_k_sample(N,pos_NV,K,n_spin,B,T2,gamma)
    x=1:N;
    x=x./N;
    S=zeros(1,length(K));
    phi=zeros(1,n_spin);
for kk=1:length(K)
        num=pos_NV(1);
        delta_x=(pos_NV(2)-pos_NV(1))/N;
        phi(1)=exp(1i*2*pi*gamma*T2*B(num));
        S(kk)=exp(1i*2*pi*K(kk)*x(num))*phi(1)*delta_x/2;
        for jj=2:n_spin-1
            num=pos_NV(jj);
            phi(jj)=exp(1i*2*pi*gamma*T2*B(num));% initial phase
            phi2=phi(jj);
            delta_x1=(pos_NV(jj)-pos_NV(jj-1))/N;
            delta_x2=(pos_NV(jj+1)-pos_NV(jj))/N;
            S(kk)=exp(1i*2*pi*K(kk)*x(num))*phi2*(delta_x1+delta_x2)/2+S(kk);
        end
        num=pos_NV(end);
        delta_x=(pos_NV(end)-pos_NV(n_spin-1))/N;
        phi(end)=exp(1i*2*pi*gamma*T2*B(num));
        S(kk)=exp(1i*2*pi*K(kk)*x(num))*phi(end)*delta_x/2+S(kk);
end
S=S./(1/n_spin);
end