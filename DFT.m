function [G1]=DFT(X,K,S,delta_k)

    G2=zeros(length(X),1);%define the Fourier series
    for jj=1:length(X)
        for ii=1:length(K)
            G2(jj)=exp(-1i*2*pi*K(ii)*X(jj))*S(ii)+G2(jj);
        end
    end
        G1=delta_k/(2*pi).*G2;
end


