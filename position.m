function [pos_NV]=position(n_spin,distance,N)
    
   
%% evenly diostributed case  
%     pos_NV=zeros(n_spin,1);
%     ran=zeros(n_spin,1);
%     for ii=0:1:n_spin-1 %NV's position relative to the magnetic field
%         pos_NV(ii+1)=ceil(1+distance*ii+sqrt(distance)*4*ran(ii+1));
%     end

 %% Totally randomly distributed case
%      pos_NV = sort(randperm(N,n_spin));%randomly distributed NV lattice 
     
%% Pertubatively randomly distributed case
%     pos_NV=zeros(n_spin,1);
%     ran=rand(n_spin,1); %randomly distributed NV lattice 
%     for ii=0:1:n_spin-1 %NV's position relative to the magnetic field
%         %pos_NV(ii+1)=ceil(1+distance*ii+sqrt(distance)*4*ran(ii+1));
%         pos_NV(ii+1)=ceil(1+distance*ii+5*ran(ii+1)-2.5); %aperture size 5 nm
%     end
%     if  pos_NV(end)>N
%         pos_NV(end)=N;
%     end
%     if pos_NV(1)<=0
%        pos_NV(1)=1;
%     end
%     

%% settled case
   %pos_NV=[350,379,532,547,582,618,645,810,871,935];
 % pos_NV=[110,170,183,194,208,225,226,230,254,301,308,401,426,432,470,841,893,914,966,991];
   pos_NV=[3,104,149,151,153,202,300,302,399,501,504,604,653,699,703,804,854,899,901,952];
   
  
end