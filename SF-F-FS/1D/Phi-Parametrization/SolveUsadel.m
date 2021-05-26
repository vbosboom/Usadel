function [phiP,phiCM] = SolveUsadel(GI,weights,T,gamma,H,phase,n_mats,maxit,tol,useprev)
    
    Delta_0 = BCSGap(T); %BCS Energy gap
    
    %initial guesses used in iteration method
    if useprev==false
        phiP = Delta_0*ones(GI.ntot,n_mats);
        phiCM = Delta_0*ones(GI.ntot,n_mats);
    else
        load SolPrev.mat
    end
    
    %solve Usadel equations for all Matsubara frequencies
    for ni =1:n_mats
        omega = (2*ni-1)*T;
        %solve Usadel equations using Picard iteration
        [phiP(:,ni),phiCM(:,ni)] = Picard(GI,weights,omega,gamma,H,phiP(:,ni),phiCM(:,ni),Delta_0,phase,tol,maxit);
    end
end