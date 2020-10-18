function [phase,phase2]=pickout(locs,PHA,reso,N)
    X_r=locs;
    phase=zeros(length(X_r),1);
    for ii=1:length(X_r)
        phase(ii)=PHA(round(X_r(ii)/reso)+1);
    end
    phase2=zeros(N,1);
    for ii=1:length(X_r)
        phase2(ceil(X_r(ii)/reso)+1)=PHA(round(X_r(ii)/reso)+1);
    end
end