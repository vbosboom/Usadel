%Solves the nonlinear matrix equation from the 1D Usadel equations in 
%the phi-parametrization using the Picard iteration
function phi = Picard(GI,weights,omega,gamma,phi,Delta_0,phase,tol,maxit)
    %iteration parameters
    Diff = 10;
    iter = 0;
    %iterate until convergence
    while (Diff>tol && iter< maxit)
        phi_old = phi; %previous iteration
        
        %build the element matrices and vectors using linear basis functions
        [S1,S2,S3,f] = BuildMatricesandVectors(GI,weights,gamma,omega,phase,phi,Delta_0);
        Stot = S1+S2+S3;
        phi = Stot\f; %update phi using the Picard iteration
        Diff = norm(phi-phi_old); %convergence criterion
        iter = iter+1;
    end
end