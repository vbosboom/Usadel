function [f1elem,f2elem,f3elem,f4elem] = GenerateVecElementVector(GI,weights,E,Theta,Chi,Gammainv,ChiGrad,h,Delta_0,phase,ind1)
    
    ChiS = [-phase/2*ones(GI.nL-1,1);zeros(GI.nS-1,1);phase/2*ones(GI.nL-1,1)];
    ThetaS = [atan(Delta_0/(-1i*E))*ones(GI.nL-1,1);zeros(GI.nS-1,1);atan(Delta_0/(-1i*E))*ones(GI.nL-1,1)];
    
    f1elem = -h/2.*(sin(Theta).*cos(Theta).*ChiGrad.*ChiGrad.*GI.PhiBS(ind1,:))*weights;
    f2elem = h/2.*(1i*E*sin(Theta).*GI.PhiBS(ind1,:))*weights;
    f3elem = h/2.*(Gammainv.*(cos(Theta).*sin(ThetaS).*cos(ChiS-Chi)-sin(Theta).*cos(ThetaS)).*GI.PhiBS(ind1,:))*weights;
    f4elem = h/2.*(Gammainv.*(sin(Theta).*sin(ThetaS).*sin(ChiS-Chi)).*GI.PhiBS(ind1,:))*weights;
end