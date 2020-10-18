%=========================
%this program is used for spectrum rid of characteristic peaks
%=========================

%% parameter
a0=5*10^(-8);%NV distance
aa0=a0*10^(6);%convert the average distance of NV into um
HW=0.8*10^(-8);%half width of the magnetic field distribution
%lamda=100;%distance (pixel size) between two peaks
B0=0.2;%amplitude of the magnetic field Gauss ,therefore phi_s(x)=2pi*gamma*T2*B_max won't exceeds 2pi, avoids the dillemma
gamma=2.803;%gyromagnetic ratio
int=10^(-9);%pixel size 1 nm
Size=10^(-6);%sample size 1000nm
n_read=Size/a0;%readout 11 NV at once
N=1000;%number of points
%Bg=20;%max gradient
pixel=int*10^6;
T2=0.3;%T2* 0.3us
C=1;%contrast 1
distance=ceil(a0/int);
n_spin=ceil(N/distance);%number of spins
space=n_spin/n_read;
delta_x=a0*10^6;%distance between NV centers in scale of um
delta_k=1/Size/10^(6);%field of view
num_pks=91;%settled due to lamda,the actual frequency would be 1000/ceil(1000/77)=76.9
lamda=ceil(N/num_pks);%distance (pixel size) between two peaks
count=50; %photon counts

%% %NV's position relative to the magnetic field
pos_NV=position(n_spin,distance,N);


%% scanning part
spacing=80; %real space scanning 50 nm
N_p=20;%sampling point
Kmax=200;%60:10:1000;
Cont=zeros(length(N_p),length(Kmax)); % used for collecting scanning data
scan_step=(N_p-1).*(spacing);

[BB]=Mag_sensor(N,B0,HW,pos_NV,int,n_spin,lamda,num_pks,0); %include scanning

for jj=1:length(N_p)
delta_K=1/spacing/N_p(jj)*1000;
start=0:spacing:spacing*N_p(jj);%moving NV sample in the real space by amount 
Recon=zeros(N_p(jj),n_spin);
scanning=zeros(N_p(jj),n_spin);


for kk=1:length(Kmax)
    
    K=0:1:Kmax(kk);

% compressed sensing(randomly tossing points)
Toss_num=Kmax-120;%ceil(Kmax(kk)/6);
Toss=sort(randperm(ceil(Kmax(kk)),Toss_num));
for ll=length(Toss):-1:1
    K(Toss(ll))=[];
end

for ii=1:N_p(jj)% real space scanning uniformly

%% plot the field and sensor
% [B]=Mag_sensor(N,B0,HW,pos_NV,int,n_spin,lamda,num_pks);
%[B]=Mag_sensor(N,B0,HW,pos_NV,int,n_spin,lamda,num_pks,start(ii)); %include scanning
B=BB(N/2+1+start(ii)+4900:3*N/2+start(ii)+4900);
%B=BB(N/2*Npixel+1+start(ii)*Npixel+4900*Npixel:3*N/2*Npixel+start(ii)*Npixel+4900*Npixel);
%% Ramsey sequence


%%%% differeny K series options
% K=sort(Kmax*rand(1,Kmax));
% K=sort(K);

% S=ksample(N,pos_NV,K,space,n_read,n_spin,B,T2,gamma);
%% add noise to S
 S=ksample_noise(N,pos_NV,K,space,n_read,n_spin,B,T2,gamma,count);

%% DFT matrix
% FFT = dftmtx(N);
% A_FFT=FFT(1:Kmax(ii)+1,:);

%% l2 norm 2020/7/21
% [G1]=L2_norm(A_FFT,S');
% ABS_norm2=abs(G1);
% ABS_norm2=ABS_norm2./max(ABS_norm2);
% PHA_norm2=angle(G1);


%% L1 norm 2020/7/22
% N = size(A_FFT,2);  %Signal dimension
% 
% cvx_begin 
% variable s(N)
% minimize( norm(A_FFT*s-S',1) ) 
% cvx_end
% 
% G1_norm=s;
% ABS_norm1=abs(G1_norm);
% ABS_norm1=ABS_norm1./max(ABS_norm1);
% PHA_norm1=angle(G1_norm);

%% compare the L2 norm and L1 norm
%     figure
%     subplot(2,1,1)
%     plot(X,ABS_norm2,'o-')
% %     hold on
% %     plot(X,ABS_norm1,'.-')
%     subplot(2,1,2)
%     plot(X,PHA,'o-')
% %     hold on
% %     plot(X,PHA_norm1,'.-')

%% perform Discrete Fourier Transformation  
    X=pixel:pixel:1;%set the lattice
    G1=DFT(X,K,S,delta_k);
    ABS=abs(G1);
    ABS=ABS./max(ABS);
    PHA=angle(G1);
    Re=real(G1);
    Im=imag(G1);
    
%% track the phase distribution after locating the NV distribution
NNN=1;
[X_r,phase]=locating(ABS,n_spin,NNN,PHA);

%% plotting
%     phi=B*(2*pi*gamma*T2);
%     ploting(N,n_spin,pos_NV,ABS,phase,phi,X,X_r,[],PHA);

%% collect scanning points
Recon(ii,1:size(phase))=phase;
X_r=X_r+start(ii)/N;
scanning(ii,1:length(X_r))=X_r;

end

%% collect real space scanning data 
YY=Recon(:);
XX=scanning(:);
XX(find(XX==0))=[];%take out the zero 
YY(find(YY==0))=[];
[XX,index]=unique(XX);
yy=zeros(length(XX),1);
for ii=1:length(XX)
    yy(ii)=YY(index(ii));
end
%[B,phi,BB]=Mag_sensor(N,B0,HW,pos_NV,int,n_spin,lamda,num_pks,start(1));
%phise=BB*(2*pi*gamma*T2);
% figure
% plot(XX,yy,'.');
% hold on
% plot(1/N:1/N:2,phise')
% xlabel('x/um')
% ylabel('amp')

%% plomb for non-uniform sampling 2020/7/8
tn=XX;
n=yy-mean(yy);

% for ii=length(tn):-1:1
%     if abs(n(ii))>0.4
%         tn(ii)=[];
%         n(ii)=[];
%     end
% end

figure
subplot(2,1,1)
plot(tn,n,'*-')
xlabel('x/um');
ylabel('phi');
title('reconstructed phase')
ofac = 2;
subplot(2,1,2)
[pxx,f]= plomb(n,tn,500,'power');

%========================this is for deleting the characteristic peaks
% Kc=1000/spacing;
% delta_f=f(2)-f(1);
% for ii=1:14
%     Pe=ceil(Kc*ii/delta_f);
%     pxx(Pe-5:Pe+5)=0;
% end

plot(f,pxx);
xlabel('K (um^{-1})');
ylabel('Power/wavevector (dB/um^{-1})')

%%
%========================this is for fourier transformation
% fs=20;
% tttt=1:20;
% gg=fft(n,20,1);
% m=abs(gg);
% f=(tttt-1)./20*fs;
% figure
% plot(f,m)
% xlabel('K (um^{-1})');
% ylabel('amp');
%% the goodness of the peak
delta_f=f(2)-f(1);
fre=N/ceil(N/num_pks);
Mid=ceil(fre/delta_f);
F_max=max(pxx(Mid-4:Mid+4));
pxx(Mid-4:Mid+4)=0;
Kc=1000/spacing;
for ii=1:30
    Pe=ceil(Kc*ii/delta_f);
    pxx(Pe-4:Pe+4)=0;
end

A_max=max(pxx);
Cont(jj,kk)=F_max/A_max;

end
end
%% 2D plot of goodness
% figure
%  plot(Kmax,Cont(1,:),'*-')
%  hold on
% % % % % plot(Kmax,Cont(2,:),'o')
% % % % % % hold on
% % % % % % plot(Kmax,Cont(3,:),'.-')
% xlabel('K_{max}');
% ylabel('goodness of the spectrum')
% grid off

% surf(scan_step,Kmax,Cont');
% grid off
% view(2)
% xlabel('overall scanning distance (nm)');
% ylabel('K_{max}');

