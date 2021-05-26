function [phiP,phiCM] = Picard(GI,weights,omega,gamma,H,phiP,phiCM,Delta_0,phase,tol,maxit)
    
    %iteration parameters
    Diff = 10;
    iter = 0;
    while (Diff>tol && iter<maxit)
        phiP_old = phiP;
        phiCM_old = phiCM;
        [S1m,S1p,S2m,S2p,S3m,S3p,fm,fp] = BuildMatricesandVectors(GI,weights,omega,gamma,H,phase,phiP,phiCM,Delta_0);
        Stotm = S1m+S2m+S3m;
        Stotp = S1p+S2p+S3p;

        phiP = Stotp\fp;
        phiCM = Stotm\fm;
        
        
        Diff = max(norm(phiCM-phiCM_old),norm(phiP-phiP_old));
        iter = iter+1;
    end
end