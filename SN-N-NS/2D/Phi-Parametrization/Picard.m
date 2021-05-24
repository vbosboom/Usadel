%Solve the Usadel equations using a Picard iteration
function [phis1,phis2,phin] = Picard(GI,points,weights,points1D,weights1D,gamma_B,gamma,omega,ksi,phis1,deltas1,phis2,deltas2,phin,phase,Delta_0,tol,maxit)
    Diff = 10;
    iter = 0;
    while ((Diff>tol) && (iter<maxit))
        phis1_old = phis1;
        phis2_old = phis2;
        phin_old = phin;
        [S1,M1,f1,S2,M2,f2,S3,M3,f3] = BuildVecMatricesandVectors(GI,points,weights,points1D,weights1D,gamma_B,gamma,omega,ksi,phis1,deltas1,phis2,deltas2,phin);
        S1tot = S1+M1;
        S2tot = S2+M2;
        S3tot = S3+M3;
        [S1tot,f1,S2tot,f2] = ProcessEssentialBoundaryConditions(GI,S1tot,f1,S2tot,f2,phase,Delta_0);
        phis1 = lsqminnorm(S1tot,f1);
        phis2 = lsqminnorm(S2tot,f2);
        phin = lsqminnorm(S3tot,f3);
        Diff = max([norm(phis1-phis1_old),norm(phis2-phis2_old),norm(phin-phin_old)]);
        iter = iter+1;
    end
end