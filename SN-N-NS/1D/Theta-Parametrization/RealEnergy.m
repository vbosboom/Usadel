%Procedure for solving the Usadel equations to find the Green's functions
function [chi,theta,success] = RealEnergy(GI,points,weights,T,phase,E,gamma,maxit,tol,useprev)
    Delta_0 = BCSGap(T); %BCS Energy gap
    if useprev==true
        load solprev.mat
    else
        %initial guesses used in iteration method
        theta = atan(Delta_0/(-1i*E))*ones(GI.ntot,1);
        chi = zeros(GI.ntot,1);
    end
    %solve the Usadel equations using the Newton method
    [chi,theta,success] = Newton(GI,points,weights,phase,E,gamma,chi,theta,Delta_0,maxit,tol);
end