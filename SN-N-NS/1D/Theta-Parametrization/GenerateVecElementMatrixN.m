function [S1elem,S2elem,JS21elem,Jf11elem,Jf12elem,Jf21elem,Jf31elem,Jf32elem,Jf41elem,Jf42elem]...
                = GenerateVecElementMatrixN(GI,weights,E,Theta,Chi,Gammainv,GradPhi,ChiGrad,h,Delta_0,phase,ind1,ind2)
    
    ChiS = [-phase/2*ones(GI.nL-1,1);zeros(GI.nS-1,1);phase/2*ones(GI.nL-1,1)];
    ThetaS = [atan(Delta_0/(-1i*E))*ones(GI.nL-1,1);zeros(GI.nS-1,1);atan(Delta_0/(-1i*E))*ones(GI.nL-1,1)];
            
    S1elem = h/2.*(GradPhi(:,ind1).*GradPhi(:,ind2))*2;
    S2elem = h/2.*(sin(Theta).^2.*GradPhi(:,ind1).*GradPhi(:,ind2))*weights;
    JS21elem = h/2.*(2*sin(Theta).*cos(Theta).*GradPhi(:,ind1).*ChiGrad.*GI.PhiBS(ind2,:))*weights;
    Jf11elem = -h/2.*((cos(Theta).^2-sin(Theta).^2).*ChiGrad.*ChiGrad.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    Jf12elem = -h/2.*(2*sin(Theta).*cos(Theta).*ChiGrad.*GI.PhiBS(ind1,:).*GradPhi(:,ind2))*weights;
    Jf21elem = h/2.*(1i*E*cos(Theta).*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    Jf31elem = h/2.*(Gammainv.*(-sin(Theta).*sin(ThetaS).*cos(ChiS-Chi)-cos(Theta).*cos(ThetaS)).*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    Jf32elem = h/2.*(Gammainv.*(cos(Theta).*sin(ThetaS).*sin(ChiS-Chi)).*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    Jf41elem = h/2.*(Gammainv.*(cos(Theta).*sin(ThetaS).*sin(ChiS-Chi)).*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    Jf42elem = -h/2.*(Gammainv.*(sin(Theta).*sin(ThetaS).*cos(ChiS-Chi)).*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
end