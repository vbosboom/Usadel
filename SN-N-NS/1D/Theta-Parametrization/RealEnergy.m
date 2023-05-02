%Procedure for solving the 1D Usadel equations in the theta,chi-parametrization
% for real energies to find the Green's functions in the normal metal
function [chi,theta,success] = RealEnergy(GI,weights,T,phase,E,gamma,maxit,tol,useprev)
    Delta_0 = BCSGap(T); %BCS Energy gap

    %load initial guess from external file
    if useprev==true
        load solprev.mat
    else
        %initial guesses used in iteration method
        theta = atan(Delta_0/(-1i*E))*ones(GI.ntot,1);
        chi = zeros(GI.ntot,1);
    end
    %solve the Usadel equations using the Newton method
    [chi,theta,success] = Newton(GI,weights,phase,E,gamma,chi,theta,Delta_0,maxit,tol);
end