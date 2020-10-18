function [AA]=recontruction(Abs_S,ABS1,N,N2,X1,X2,phi)
     S=zeros(1,500);
     S(1)=(Abs_S(1)*(X1(1)/X2(1)))/2/N2;
     S(26:36)=Abs_S(2:12)./N2*2;
     AA=zeros(N,1);
    AA(1:N/2)=S;
    AA(N/2+1)=0;
    AA(N/2+2:N+1)=fliplr(S);
    AA(N+1)=[];
     


     Y1=fft(ABS1,N);
     G1 = (Y1/N*2);%amp need to be divided by N/2
     G1(1)=G1(1)/2;% zero mode point need to be divided by N
     G1=abs(G1);
     f1=(1/N:1/N:1);

     Y2=fft(AA,N);
     G2 = (Y2/N*2);%amp need to be divided by N/2
     G2(1)=G2(1)/2;% zero mode point need to be divided by N
    G2=abs(G2);

figure
plot(f1,G1*max(phi)/max(G1))
hold on
plot(f1,G2*max(phi)/max(G1));
title('reconstructed field');
xlabel('x/um');
ylabel('amp')
 legend('input field','reconstructed field')    
    
end