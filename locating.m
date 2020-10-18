function [X_r,phase]=locating(ABS,n_spin,NNN,PHA)
%      [pks,locs] = findpeaks(ABS);
%     for ii=length(locs):-1:1
%         if pks(ii)<=0.35
%             locs(ii)=[];
%         end
%     end
% locs=[3,104,149,151,153,202,300,301,399,502,504,604,653,698,704,804,854,899,902,952];
 locs=[3,104,202,399,604,653,804,854,952];
     X_r=locs;

    phase=zeros(length(X_r),1);
    for ii=1:length(X_r)
        phase(ii)=PHA(X_r(ii));
    end
    X_r=locs./1000;
    %X_r=locs./n_spin/NNN-1/n_spin/NNN;% convert the scale to um
end