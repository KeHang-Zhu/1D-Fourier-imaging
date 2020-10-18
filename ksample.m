function [S]=ksample(N,pos_NV,K,space,n_read,n_spin,B,T2,gamma)%N is the multiple of 
    x=1:N;
    x=x./1000;%change into um
S=zeros(1,length(K));
    phi=zeros(n_spin,1);
for kk=1:length(K)
        for jj=1:n_spin %readout 11 sensor at once
            num=pos_NV(jj);
            %              phi((ii-1)*n_read+jj)=exp(1i*2*pi*gamma*T2*B(num));
            phi(jj)=exp(1i*2*pi*gamma*T2*B(num));% initial phase
            phi2=phi(jj);
            S(kk)=exp(1i*2*pi*K(kk)*x(num))*phi2+S(kk);
            %              S(ii,kk)=exp(1i*2*pi*(K(kk)*x(num)))*B(num)+S(ii,kk);
        end
end

end