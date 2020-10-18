function [B1,phi]=Mag_sensor(N,B0,HW,pos_NV,int,n_spin,lamda,num_pks,start)


    x=(1:N);
    x=x/1000;%change into um
    X0=ceil(N/num_pks/2);
    
%     for ii=1:N
%         %       B(ii)=B0/(1+((ii-ceil(2*N/3))*int/HW)^2)+B0/(1+((ii-ceil(N/3))*int/HW)^2);%two lorentzian function distribution
%         B(ii)=B0/(1+((ii-N/2)*int/HW)^2)+B0/(1+((ii-N/2+200)*int/HW)^2);%+B0/(1+((ii-N/2+100))*int/HW)^2;
%     end


    B1=zeros(10*N,1);
    R=zeros(10*N,1);
for jj=1:num_pks*10
     X=X0+(jj-1)*lamda;
    for ii=1:N*10
        R(ii)=B0/(1+((ii-X)*int/HW)^2);
        % R(ii)=B0/(1+((ii-X)*int/HW)^2);%+B0/(1+((ii-N/2+100))*int/HW)^2;%Lorentzian
%          R(ii)=B0*exp(-((ii-X)*int/HW)^2);%Gaussian
    end
    B1=B1+R;
end

B1=B1-mean(B1);
B1=B1.*3+0.2;

% B=B1(N/2+1+start:3*N/2+start)-mean(B1(N/2+1+start:3*N/2+start));
% B=B.*3+0.2;
% BB=B1(N/2+1+start:5*N/2+start)-mean(B1(N/2+1+start:5*N/2+start));
% BB=BB.*3+0.2;
phi=2*pi*0.3*2.803.*B1;


%% plot single peak

    jj=num_pks;
    X=X0+(jj-1)*lamda;
    ii=1:2*N;
     B2=B0./(1+((ii-X).*int./HW).^2);%lorentzian
%      B2=B0.*exp(-((ii-X).*int./HW).^2);%Gaussian
    
Bsingle=B2(N/2+1:3*N/2);
phi_single=2*pi*0.3*2.803.*Bsingle;

%% plot adjancent peak
    jj=num_pks+1;
    X=X0+(jj-1)*lamda;
    ii=1:2*N;
     B2=B0./(1+((ii-X).*int./HW).^2);%lorentzian
%     B2=B0*exp(-((ii-X).*int./HW).^2);%Gaussian

Bsingle2=B2(N/2+1:3*N/2);
phi_single2=2*pi*0.3*2.803.*Bsingle2;

%     B=zeros(N,1);
%     X0=floor(N/num_pks/2);
%     n=floor(N/num_pks);
% for ii=1:num_pks
%     for jj=1:n
%         B((ii-1)*n+jj)=B0/(1+((jj-X0)*int/HW)^2);
%     end
% end

%% sin field
% ii=1:N;
% B=B0*sin(2*pi*ii/N*33)+B0;
% phi=2*pi*0.3*2.803.*B;
%% triangular field
% B=zeros(N,1);
% for ii=1:num_pks
%     for jj=1:round(lamda/2)
%         B((ii-1)*lamda+jj)=B0*(jj/lamda*2)+B0;
%     end 
%     for jj=round(lamda/2)+1:lamda
%         B((ii-1)*lamda+jj)=-B0*(jj/lamda*2)+3*B0;
%     end
% end
% B(N+1:end)=[];
% phi=2*pi*0.3*2.803.*B;
%% periodic domain wall
% B=zeros(N,1);
% for ii=1:num_pks
%     for jj=1:round(lamda/2)
%         B((ii-1)*lamda+jj)=B0;
%     end 
%     for jj=round(lamda/2)+1:lamda
%         B((ii-1)*lamda+jj)=-B0;
%     end
% end
% B(N+1:end)=[];
% phi=2*pi*0.3*2.803.*B;


%% plotting

% figure
% plot(x,phi)
% hold on
% plot(x,phi_single,'b.-');
% hold on
% plot(x,phi_single2,'r.-');
% xlabel('x/um');
% ylabel('phi');
% title('Field ditribution and sensor location');
% 
% hold on
% for ii=1:1:n_spin
%     plot(x(pos_NV(ii)),0,'r*')
% end
% hold off
% legend('multi-peaks','single-peak','sensor locs')
%     
%%    
%     figure
%     plot(x,phi)
%     hold on
%     plot(x,phi_single,'.-');
%     xlabel('x/um');
%     ylabel('B/Gauss');
%     title('Field ditribution and sensor location');
% 
%     hold on
%     for ii=1:1:n_spin
%         plot(x(pos_NV(ii)),0,'r*')
%     end
%     hold off
    
%     figure
%     plot(x,B);
%     hold on
%     plot(x,B2,'--');
%     xlabel('x/um');
%     ylabel('B/Gauss');
%     title('Field ditribution and sensor location');
%     hold on
%     
%     for ii=1:1:n_spin
%         plot(x(pos_NV(ii)),0,'r*')
%     end
%     hold off


end