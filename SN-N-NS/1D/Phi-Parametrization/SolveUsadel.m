%Procedure for solving the 1D-Usadel equations in the phi-parametrization
%to find the Green's functions in the normal metal layer
function phin = SolveUsadel(GI,weights,T,gamma,phase,n_mats,maxit,tol)
    Delta_0 = BCSGap(T); %BCS Energy gap
    
    %initial guesses used in iteration method
    phin = zeros(GI.ntot,n_mats); %Container for the Green's function solutions at different Matsubara frequencies

    %solve Usadel equations for all Matsubara frequencies
    for ni = 1:n_mats
        omega = (2*ni-1)*T; %current matsubara frequency
        
        %solve Usadel equations using the Picard method
        phin(:,ni) = Picard(GI,weights,omega,gamma,phin(:,ni),Delta_0,phase,tol,maxit);
    end
end