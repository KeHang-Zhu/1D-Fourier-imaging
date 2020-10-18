function [S]=ksample_noise(N,pos_NV,K,space,n_read,n_spin,B,T2,gamma,max)%N is the multiple of 
    x=1:N;
    x=x./N;
S=zeros(1,length(K));

    %phi=zeros(n_spin,1);
   
for kk=1:length(K)
    for ii=1:space %along the line
        for jj=1:n_read %readout 11 sensor at once
            num=pos_NV((ii-1)*n_read+jj);
           
            %phi((ii-1)*n_read+jj)=exp(1i*2*pi*gamma*T2*B(num));% initial phase
            %phi2=phi((ii-1)*n_read+jj);
            %S(ii,kk)=exp(1i*2*pi*K(kk)*x(num))*phi2+S(ii,kk);
           
            
            aa=max*cos(2*pi*K(kk)*x(num)+2*pi*gamma*T2*B(num));
            aa=floor(aa);
            bb=max*sin(2*pi*K(kk)*x(num)+2*pi*gamma*T2*B(num));
            bb=floor(bb);
            if aa>0
                Sx=poissrnd(aa)/max;
            elseif aa<0
                Sx=-poissrnd(-aa)/max;
            end
            if bb>0
                Sy=poissrnd(bb)/max;
            elseif bb<0
                Sy=-poissrnd(-bb)/max;
            end  
            S(ii,kk)=S(ii,kk)+Sx+1i*Sy;
        end
    end
end
  

end