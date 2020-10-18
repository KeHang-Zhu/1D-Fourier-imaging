function []=ploting(N,n_spin,pos_NV,ABS,phase,phi,X,X_r,phase2,PHA)
    x=1:N;
    x=x./1000;%convert into um
    pos=zeros(N,1);
    for ii=1:1:n_spin
        pos(pos_NV(ii))=1;
    end
    figure
    subplot(2,1,1)
    plot(x,pos,'*-');
    hold on
    title('NV distribution');
    xlabel('x/um');
    plot(X,ABS,'o-');


    subplot(2,1,2)
    plot(x,phi)
    hold on
     plot(X_r,phase,'*')
%     plot(X_r,phase2,'*')
%      plot(pos_NV./1000,phase2,'o')
%     plot(pos_NV./1000,phase2,'o')
       plot(X,PHA,'.-');
    xlabel('x/um');
    ylabel('phi');
    title('phi distribution');
    legend('input field','reconstructed','F^{-1}(S)');
    
%     X_actual=pos_NV./N;
%     Y=X_r-X_actual;
%     pp=1:n_spin;
%     figure
%     plot(pp,Y);
%     xlabel('NV number')
%     ylabel('difference between actual and extracted postion')
end