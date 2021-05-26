function [gammaF,gammaTF,succes] = Newton(GI,weights,omega,gammaBM,H,gammaF,gammaTF,Delta_0,phase,tol,maxit)
    
    succes = true;
    %iteration parameters
    Diff = 10;
    iter = 0;
    while (Diff>tol && iter<maxit)
        [S1,S1T,S2,S2T,S3,S3T,S4,S4T,S5,S5T,S6,S6T,f,fT,...
        JS2NN,JS2NT,JS2TN,JS2TT,JS6NN,JS6TT] = BuildMatricesandVectorsN(GI,weights,omega,gammaBM,H,phase,gammaF,gammaTF,Delta_0);
        
        Stot = S1+S2+S3+S4+S5+S6;
        StotT = S1T+S2T+S3T+S4T+S5T+S6T;
        
        h = [Stot*gammaF-f;StotT*gammaTF-fT];
        Jac = [Stot+JS2NN+JS6NN,JS2NT;JS2TN,StotT+JS2TT+JS6TT];
        
        %apply linesearch to the Newton method
        pk = lsqminnorm(Jac,h);
        alpha = linesearch(0.5,10,pk,norm(h),gammaF,gammaTF,GI,weights,omega,gammaBM,H,phase,Delta_0);
        gammaF = gammaF-alpha*pk(1:GI.ntot);
        gammaTF = gammaTF-alpha*pk(GI.ntot+1:end);
        
        Diff = norm(h);
        iter = iter+1;
    end
    if iter>=maxit
        warning('Newton method did not converge')
        succes = false;
    end
end