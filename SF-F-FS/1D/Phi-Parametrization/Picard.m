function [phiP,phiCM] = Picard(GI,weights,omega,gamma,H,phiP,phiCM,Delta_0,phase,tol,maxit)
    
    %Iteration parameters
    Diff = 10; %Current difference between solutions
    iter = 0; %Iteration counter

    while (Diff>tol && iter<maxit)
        %Save previous iteration
        phiP_old = phiP;
        phiCM_old = phiCM;
        
        %Assemble the finite element matrices and vectors
        [S1m,S1p,S2m,S2p,S3m,S3p,fm,fp] = BuildMatricesandVectors(GI,weights,omega,gamma,H,phase,phiP,phiCM,Delta_0);
        Stotm = S1m+S2m+S3m;
        Stotp = S1p+S2p+S3p;

        %Apply the Picard iteration
        phiP = Stotp\fp;
        phiCM = Stotm\fm;
        
        %Calculate norm of differences
        Diff = max(norm(phiCM-phiCM_old),norm(phiP-phiP_old));
        iter = iter+1;
    end
end