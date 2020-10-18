% function [G1,ABS2,SSS,fff1,fff2,fit_info]=comparison(phi,phase,Abs_S,K,distance,N,NNN)
function [G1,ABS2,SSS,Ks]=comparison(phi,phase,Abs_S,K,distance,N,NNN,delta_k)
Y1=fft(phi,N);
G1 = (Y1/N*2);%amp need to be divided by N/2
 G1(1)=G1(1)/2;% zero mode point need to be divided by N
f1=(0:delta_k:(N-1)*delta_k);
   % G2=abs(Y3/point);
   ABS1=abs(G1);
%     IMM1=real(G1);
%     REL1=imag(G1);


N2=round(N/distance);
Y2=fft(phase,N2);
G2 = (Y2/N2*2);%amp need to be divided by N/2
 G2(1)=G2(1)/2;% zero mode point need to be divided by N
f2=(0:1:N2-1);
% G2=abs(Y3/point);
ABS2=abs(G2);
%     IMM2=real(G2);
%     REL2=imag(G2);

%% curve fitting
% FF1=ABS1(2:N/2);
% FF1=FF1';
% XX1=f1(2:N/2);
% fun = @(x,xdata)x(1)./(1+(xdata-x(2)).^2./x(3)^2);
% %fun = @(x,xdata)x(1).*exp(-(xdata-x(2)).^2./x(3)^2);
%     x0 = [0.05,34,1.14];%fitting initial value
%     fff1 = lsqcurvefit(fun,x0,XX1,FF1);
%     X=0:0.01:f1(N/2);
%     fphi=fff1(1)./(1+(X-fff1(2)).^2./fff1(3)^2);
%     
% FF=Abs_S(2:floor(length(K)/2/NNN))./N2;
% XX=K(2:floor(length(K)/2/NNN));
% 
%     [pks,locs]=findpeaks(FF);
%     x1 = [pks,locs,1];%fitting initial value
%       fff2 = lsqcurvefit(fun,x1,XX,FF);
%     x2=0:0.01:XX(end);
%     fphi_c=fff2(1)./(1+(x2-fff2(2)).^2./fff2(3)^2);
%     fit_info=[fff1,fff2]
    
    
Abs_S(2:end)=Abs_S(2:end)./N2.*2;%amp need to be divided by N/2
Abs_S(1)=Abs_S(1)./N2;    % zero mode point need to be divided by N
figure
plot(f1(1:N/2),ABS1(1:N/2),'*');

SSS=Abs_S(1:round(length(K)/NNN/2));
% hold on
% plot(f2(1:N2/2),ABS2(1:N2/2),'o');
% hold on
% plot(K(1:round(length(K)/NNN/2)),SSS,'+-');

%% revised 2020/7/17
% hold on
% plot(f2,ABS2,'o');
% hold on
% plot(K,Abs_S,'+-');
%%


% hold on
% plot(X,fphi,'r-');
% hold on
% plot(x2,fphi_c,'b-');
xlabel('K/um^{-1}');
ylabel('Amp');
 legend('F[phi]','F[phi_{meas}]','S(K)');
% legend('F[phi]','F[phi_{meas}]','S(K)');%,'fitting','fitting');
%  legend('F[phi]','fitting','S(K)','fitting');

Ks=K(round(length(K)/NNN/2));

end