%solves the nonlinear system from FEM using the Picard iteration
function phi = Picard(GI,weights,omega,gamma,phi,Delta_0,phase,tol,maxit)
    %iteration parameters
    Diff = 10;
    iter = 0;
    %iterate untill convergence
    while (Diff>tol && iter< maxit)
        phi_old = phi; %previous iteration
        
        %build the element matrices and vectors in FEM
        [S1,S2,S3,f] = BuildMatricesandVectors(GI,weights,gamma,omega,phase,phi,Delta_0);
        Stot = S1+S2+S3;
        phi = Stot\f;
        Diff = norm(phi-phi_old);
        iter = iter+1;
    end
end