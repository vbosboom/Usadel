function [deltas1,deltas2,fs1,gs1,Fs1,Gs1,fs2,gs2,Fs2,Gs2] = DeltaIteration(deltas1,deltas2,Sumfs1,Sumfs2,SumOmega,T,iter2,fs1,gs1,Fs1,Gs1,fs2,gs2,Fs2,Gs2,M)
    if iter2==0
        gs1(:,end) = (T*Sumfs1)./(log(T)+2*T*SumOmega);
        fs1(:,end) = gs1(:,end)-deltas1;
        deltas1 = gs1(:,end); %simple Picard iteration
        
        gs2(:,end) = (T*Sumfs2)./(log(T)+2*T*SumOmega);
        fs2(:,end) = gs2(:,end)-deltas2;
        deltas2 = gs2(:,end); %simple Picard iteration
        return
    end
    
    c = min(iter2,M);
    fs1 = circshift(fs1,-1,2);
    gs1 = circshift(gs1,-1,2);
    Fs1 = circshift(Fs1,-1,2);
    Gs1 = circshift(Gs1,-1,2);
    
    fs2 = circshift(fs2,-1,2);
    gs2 = circshift(gs2,-1,2);
    Fs2 = circshift(Fs2,-1,2);
    Gs2 = circshift(Gs2,-1,2);
    
    gs1(:,end) = (T*Sumfs1)./(log(T)+2*T*SumOmega);
    fs1(:,end) = gs1(:,end)-deltas1;
    Fs1(:,end) = fs1(:,end)-fs1(:,end-1);
    Gs1(:,end) = gs1(:,end)-gs1(:,end-1);
    
    gs2(:,end) = (T*Sumfs2)./(log(T)+2*T*SumOmega);
    fs2(:,end) = gs2(:,end)-deltas2;
    Fs2(:,end) = fs2(:,end)-fs2(:,end-1);
    Gs2(:,end) = gs2(:,end)-gs2(:,end-1);
    
    gammas1 = lsqlin(Fs1(:,end-c+1:end),fs1(:,end));
    gammas2 = lsqlin(Fs2(:,end-c+1:end),fs2(:,end));
    deltas1 = gs1(:,end)-Gs1(:,end-c+1:end)*gammas1;
    deltas2 = gs2(:,end)-Gs2(:,end-c+1:end)*gammas2;
end