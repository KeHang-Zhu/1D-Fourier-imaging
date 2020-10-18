HW=1*10^(-2);
K=0:0.1:N/2;
delta_x=lamda/N;
f=zeros(length(K),1);
f(2:end)=T2*gamma*B0*pi*HW*abs((1-exp(-1i.*2*pi*K(2:end)*(num_pks+1)*delta_x))./(1-exp(-1i.*2*pi*K(2:end)*delta_x))).*exp(-HW^2.*pi^2*K(2:end).^2);
f(1)=T2*gamma*B0*pi*HW*(num_pks+1);
figure
plot(K,f);
xlabel('K/um^{-1}');
ylabel('Amp');
title('Fourier transformation of a periodic Lorentzian function');
