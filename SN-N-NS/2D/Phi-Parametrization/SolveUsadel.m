function [phis1tot,phis2tot,phintot,deltas1,deltas2,Difference] = SolveUsadel(GI,points,weights,points1D,weights1D,T,gamma_B,gamma,phase,ksi,n_mats,maxit,tol,useprev)
    Delta_0 = BCSGap(T); %Greens function in bulk superconductor

    if useprev==true
        load SolPrev.mat
    else
        %initial guesses for iteration procedure
        deltas1 = ones(GI.ns1,1);
        phis1 = Delta_0*exp(-1i*phase/2)*ones(GI.ns1,1);
        deltas2 = ones(GI.ns2,1);
        phis2 = Delta_0*exp(1i*phase/2)*ones(GI.ns2,1);
        phin = 0.125*ones(GI.nn,1);
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
    phis1tot = zeros(GI.ns1,n_mats);
    phis2tot = zeros(GI.ns2,n_mats);
    phintot = zeros(GI.nn,n_mats);
    while (Diff2>tol && iter2<maxit)
        deltas1_old = deltas1;
        deltas2_old = deltas2;
        Sumfs1 = zeros(GI.ns1,1);
        Sumfs2 = zeros(GI.ns2,1);
        SumOmega = 0;
        for ni=1:n_mats
            omega = (2*ni-1)*T;
            [phis1,phis2,phin] = Picard(GI,points,weights,points1D,weights1D,gamma_B,gamma,omega,ksi,phis1,deltas1,phis2,deltas2,phin,phase,Delta_0,tol,maxit);
            phis1tot(:,ni) = phis1;
            phis2tot(:,ni) = phis2;
            phintot(:,ni) = phin;
            SumOmega = SumOmega+1/omega;
            Sumfs1 = Sumfs1+phis1./(sqrt(omega^2+abs(phis1).^2));
            Sumfs2 = Sumfs2+phis2./(sqrt(omega^2+abs(phis2).^2));
        end
        [deltas1,deltas2,fs1,gs1,Fs1,Gs1,fs2,gs2,Fs2,Gs2] = DeltaIteration(deltas1,deltas2,Sumfs1,Sumfs2,SumOmega,T,iter2,fs1,gs1,Fs1,Gs1,fs2,gs2,Fs2,Gs2,M);
        Diff2 = max([norm(deltas1-deltas1_old),norm(deltas2-deltas2_old)]);
        iter2 = iter2+1;
        Difference(iter2) = Diff2;
    end
end