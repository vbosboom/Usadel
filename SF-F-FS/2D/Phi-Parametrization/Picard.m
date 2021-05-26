function [phis1P,phis1CM,phis2P,phis2CM,phifP,phifCM] = Picard(GI,points,weights,points1D,weights1D,gamma_B,gamma,omega,ksi,phis1P,phis1CM,deltas1,phis2P,phis2CM,deltas2,phifP,phifCM,H,phase,Delta_0,tol,maxit)
    
    Diff = 10;
    iter = 0;
    while ((Diff>tol) && (iter<maxit))
        phis1P_old = phis1P;
        phis1CM_old = phis1CM;
        phis2P_old = phis2P;
        phis2CM_old = phis2CM;
        phifP_old = phifP;
        phifCM_old = phifCM;
        [S1m,S1p,M1m,M1p,f1m,f1p,S2m,S2p,M2m,M2p,f2m,f2p,S3m,S3p,M3m,M3p,f3m,f3p] = BuildVecMatricesandVectors(GI,points,weights,points1D,weights1D,gamma_B,gamma,omega,ksi,phis1P,phis1CM,deltas1,phis2P,phis2CM,deltas2,phifP,phifCM,H);
        S1mtot = S1m+M1m;
        S1ptot = S1p+M1p;
        S2mtot = S2m+M2m;
        S2ptot = S2p+M2p;
        S3mtot = S3m+M3m;
        S3ptot = S3p+M3p;
        [S1mtot,f1m,S1ptot,f1p,S2mtot,f2m,S2ptot,f2p] = ProcessEssentialBoundaryConditions(GI,S1mtot,f1m,S1ptot,f1p,S2mtot,f2m,S2ptot,f2p,phase,Delta_0);
        phis1CM = S1mtot\f1m;
        phis1P = S1ptot\f1p;
        phis2CM = S2mtot\f2m;
        phis2P = S2ptot\f2p;
        phifCM = S3mtot\f3m;
        phifP = S3ptot\f3p;
        Diff = max([norm(phis1CM-phis1CM_old),norm(phis1P-phis1P_old),norm(phis2CM-phis2CM_old),norm(phis2P-phis2P_old),norm(phifCM-phifCM_old),norm(phifP-phifP_old)]);
        iter = iter+1;
    end
end