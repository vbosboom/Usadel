%calculates the Green's functions theta and chi in all regions for the
%given input geometry parameters and energy gaps deltas1,deltas2

function [thetas1,chis1,thetas2,chis2,thetan,chin,success] = RealEnergy(GI,points,weights,points1D,weights1D,T,gamma_B,gamma,phase,ksi,E,deltas1,deltas2,itermax,tol,useprev)
    
    Delta_0 = BCSGap(T); %BCS energy gap

    %initial guesses
    if useprev==false
        thetas1 = atan(Delta_0/(-1i*E))*ones(GI.ns1,1);
        chis1 = -phase/2*ones(GI.ns1,1);
        thetas2 = atan(Delta_0/(-1i*E))*ones(GI.ns2,1);
        chis2 = phase/2*ones(GI.ns2,1);
        thetan = atan(Delta_0/(-1i*E))*ones(GI.nn,1); %option 1
        chin = zeros(GI.nn,1);
    else   
        load solprev.mat
    end
    
    %calculate solution using the Newton method
    [thetas1,chis1,thetas2,chis2,thetan,chin,success] = Picard(GI,points,weights,points1D,weights1D,phase,E,gamma_B,gamma,ksi,thetas1,chis1,thetas2,chis2,thetan,chin,deltas1,deltas2,Delta_0,itermax,tol);
end