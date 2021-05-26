function [gammaF,gammaTF,succes] = SolveReal(GI,weights,T,gammaBM,H,phase,E,maxit,tol,useprev)
    
    Delta_0 = BCSGap(T); %BCS Energy gap
    
    %initial guesses used in iteration method
    if useprev==false
        gammaF = zeros(GI.ntot,1);
        gammaTF = -conj(gammaF);
    else
        load SolPrevReal.mat
    end
    %solve Usadel equations using Picard iteration
    [gammaF,gammaTF,succes] = Newton(GI,weights,-1i*E,gammaBM,H,gammaF,gammaTF,Delta_0,phase,tol,maxit);
end