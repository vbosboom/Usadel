%Compute the vectors obtained in the weak formulation of the 1D Usadel 
%equation in the theta,chi-parametrization
%f1 corresponds to -(sin(theta_N)*cos(theta_N)*(dchi(theta_N)d_x)^2
%f2 corresponds to iE*sin(theta_N)
%f3 corresponds to sin(theta_N)*sin(theta_S)*cos(chi_S-Chi_N)-sin(theta_N)*cos(theta_S))/gamma_{BM}
%f4 corresponds to (1/gamma_{BM})*sin(theta_N)*sin(theta_S)*sin(chi_S-chi_N)

function [f1elem,f2elem,f3elem,f4elem] = GenerateElementVector(GI,weights,E,Theta,Chi,Gammainv,ChiGrad,h,Delta_0,phase,ind1)
    
    ChiS = [-phase/2*ones(GI.nL-1,1);zeros(GI.nS-1,1);phase/2*ones(GI.nL-1,1)];
    ThetaS = [atan(Delta_0/(-1i*E))*ones(GI.nL-1,1);zeros(GI.nS-1,1);atan(Delta_0/(-1i*E))*ones(GI.nL-1,1)];
    
    f1elem = -h/2.*(sin(Theta).*cos(Theta).*ChiGrad.*ChiGrad.*GI.PhiBS(ind1,:))*weights;
    f2elem = h/2.*(1i*E*sin(Theta).*GI.PhiBS(ind1,:))*weights;
    f3elem = h/2.*(Gammainv.*(cos(Theta).*sin(ThetaS).*cos(ChiS-Chi)-sin(Theta).*cos(ThetaS)).*GI.PhiBS(ind1,:))*weights;
    f4elem = h/2.*(Gammainv.*(sin(Theta).*sin(ThetaS).*sin(ChiS-Chi)).*GI.PhiBS(ind1,:))*weights;
end