function [S1mtot,f1m,S1ptot,f1p,S2mtot,f2m,S2ptot,f2p] = ProcessEssentialBoundaryConditions(GI,S1mtot,f1m,S1ptot,f1p,S2mtot,f2m,S2ptot,f2p,phase,Delta_0)
    

    %Essential boundary conditions for S1
    for i = 1:length(GI.ess_bc_index_s1)
        k = GI.convs1(GI.ess_bc_index_s1(i));
        f1m = f1m- (S1mtot(:,k))*Delta_0*exp(1i*phase/2);
        f1p = f1p- (S1ptot(:,k))*Delta_0*exp(-1i*phase/2);
    end
    
    for i = 1:length(GI.ess_bc_index_s1)
        k = GI.convs1(GI.ess_bc_index_s1(i));
        S1mtot(:,k) = 0;
        S1mtot(k,:) = 0;
        S1mtot(k,k) = 1;
        f1m(k) = Delta_0*exp(1i*phase/2);
        S1ptot(:,k) = 0;
        S1ptot(k,:) = 0;
        S1ptot(k,k) = 1;
        f1p(k) = Delta_0*exp(-1i*phase/2);
    end
    
    %Essential boundary conditions for S2
    for i = 1:length(GI.ess_bc_index_s2)
        k = GI.convs2(GI.ess_bc_index_s2(i));
        f2m = f2m- (S2mtot(:,k))*Delta_0*exp(-1i*phase/2);
        f2p = f2p- (S2ptot(:,k))*Delta_0*exp(1i*phase/2);
    end

    for i = 1:length(GI.ess_bc_index_s2)
        k = GI.convs2(GI.ess_bc_index_s2(i));
        S2mtot(:,k) = 0;
        S2mtot(k,:) = 0;
        S2mtot(k,k) = 1;
        f2m(k) = Delta_0*exp(-1i*phase/2);
        S2ptot(:,k) = 0;
        S2ptot(k,:) = 0;
        S2ptot(k,k) = 1;
        f2p(k) = Delta_0*exp(1i*phase/2);
    end
end