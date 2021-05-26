function [phis1P,phis1CM,phis2P,phis2CM,phifP,phifCM,deltas1,deltas2,Difference] = SolveUsadel(GI,points,weights,points1D,weights1D,T,gamma_B,gamma,H,phase,ksi,n_mats,maxit,tol,useprev)

    Delta_0 = BCSGap(T); %BCS Energy gap
    
    if useprev==true
        load SolPrev.mat
    else
        %initial guesses for iteration procedure
        deltas1 = Delta_0*exp(-1i*phase/2)*ones(GI.ns1,1);
        phis1P = Delta_0*exp(-1i*phase/2)*ones(GI.ns1,n_mats);
        phis1CM = Delta_0*exp(1i*phase/2)*ones(GI.ns1,n_mats);
        deltas2 = Delta_0*exp(1i*phase/2)*ones(GI.ns2,1);
        phis2P = Delta_0*exp(1i*phase/2)*ones(GI.ns2,n_mats);
        phis2CM = Delta_0*exp(-1i*phase/2)*ones(GI.ns2,n_mats);
        phifP = 0.125*ones(GI.nn,n_mats);
        phifCM = 0.125*ones(GI.nn,n_mats);
    end
    
    %sum over all used matsubara frequencies
    Diff2 = 10;
    iter2 = 0;
    Difference = 0;
    
    %Delta anderson acceleration
    M = 3;
    fs1 = zeros(GI.ns1,M+1);
    gs1 = zeros(GI.ns1,M+1);
    Fs1 = zeros(GI.ns1,M);
    Gs1 = zeros(GI.ns1,M);
    
    fs2 = zeros(GI.ns2,M+1);
    gs2 = zeros(GI.ns2,M+1);
    Fs2 = zeros(GI.ns2,M);
    Gs2 = zeros(GI.ns2,M);
    
    while (Diff2>tol && iter2<maxit)
        deltas1_old = deltas1;
        deltas2_old = deltas2;
        Sumfs1 = zeros(GI.ns1,1);
        Sumfs2 = zeros(GI.ns2,1);
        SumOmega = 0;
        
        for ni=1:n_mats
            omega = (2*ni-1)*T;
            [phis1P(:,ni),phis1CM(:,ni),phis2P(:,ni),phis2CM(:,ni),phifP(:,ni),phifCM(:,ni)] = Picard(GI,points,weights,points1D,weights1D,gamma_B,gamma,omega,ksi,phis1P(:,ni),phis1CM(:,ni),deltas1,phis2P(:,ni),phis2CM(:,ni),deltas2,phifP(:,ni),phifCM(:,ni),H,phase,Delta_0,tol,maxit);
            SumOmega = SumOmega+1/omega;
            Sumfs1 = Sumfs1+phis1P(:,ni)./(sqrt(omega^2+phis1P(:,ni).*phis1CM(:,ni)));
            Sumfs1 = Sumfs1+conj(phis1CM(:,ni))./(sqrt(omega^2+conj(phis1CM(:,ni)).*conj(phis1P(:,ni))));
            Sumfs2 = Sumfs2+phis2P(:,ni)./(sqrt(omega^2+phis2P(:,ni).*phis2CM(:,ni)));
            Sumfs2 = Sumfs2+conj(phis2CM(:,ni))./(sqrt(omega^2+conj(phis2CM(:,ni)).*conj(phis2P(:,ni))));
        end
        [deltas1,deltas2,fs1,gs1,Fs1,Gs1,fs2,gs2,Fs2,Gs2] = DeltaIteration(deltas1,deltas2,Sumfs1,Sumfs2,SumOmega,T,iter2,fs1,gs1,Fs1,Gs1,fs2,gs2,Fs2,Gs2,M);
        Diff2 = max([norm(deltas1-deltas1_old),norm(deltas2-deltas2_old)]);
        iter2 = iter2+1;
        Difference(iter2) = Diff2;
    end
end