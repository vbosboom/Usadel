%Solves the nonlinear matrix equation from the 2D Usadel equations in 
%the phi-parametrization using the Picard iteration
function [phis1,phis2,phin] = Picard(GI,points,weights,points1D,weights1D,gamma_B,gamma,omega,ksi,phis1,deltas1,phis2,deltas2,phin,phase,Delta_0,tol,maxit)
    %Iteration parameters
    Diff = 10;
    iter = 0;
    %Iterate until covergence
    while ((Diff>tol) && (iter<maxit))
        %prevous iteration
        phis1_old = phis1;
        phis2_old = phis2;
        phin_old = phin;

        %build the element matrices and vectors using linear basis functions
        [S1,M1,f1,S2,M2,f2,S3,M3,f3] = BuildMatricesandVectors(GI,points,weights,points1D,weights1D,gamma_B,gamma,omega,ksi,phis1,deltas1,phis2,deltas2,phin);
        S1tot = S1+M1;
        S2tot = S2+M2;
        S3tot = S3+M3;

        %Apply the boundary conditions
        [S1tot,f1,S2tot,f2] = ProcessEssentialBoundaryConditions(GI,S1tot,f1,S2tot,f2,phase,Delta_0);
        
        %calculate new solution using a Picard step
        phis1 = lsqminnorm(S1tot,f1);
        phis2 = lsqminnorm(S2tot,f2);
        phin = lsqminnorm(S3tot,f3);
        
        %Convergence criterion
        Diff = max([norm(phis1-phis1_old),norm(phis2-phis2_old),norm(phin-phin_old)]);
        iter = iter+1;
    end
end